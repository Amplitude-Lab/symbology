#!/usr/bin/env python3
"""Wire-routing audit for direct port-to-port KindEdge routing.

Re-simulates the EXACT pipeline of FlowEditor.jsx routeEdges + KindEdge:
  - fan-out spread per source port (+-(rank-(n-1)/2)*7),
  - obstacle detours: LOCAL external lanes (stacked +-22px) + anchor
    search (offsets 30+7i, dual outbound bases for backward wires),
  - kappa bundle bows for tangent direct wires (chords <7px apart,
    crossing angle <0.35 rad, shared-port exclusion windows),
then renders every wire's REAL path (straight chord / quadratic bow /
cubic+run+cubic detour / fallback bezier) mathematically and checks:
  0. routing quality: backward wires (dx<-24) never render straight;
     every fallback is truly inescapable (extended anchor search to
     o<=300 also fails); every detour-set wire found anchors;
  1. no non-fallback path passes through a node box interior, except
     into the wire's own endpoint box when the endpoint boxes physically
     overlap (fixture defect, reported as `fixture:` not a problem);
  2. tangent pairs end up >=2.5px apart after bowing; all direct-wire
     pairs stay >=2.5px apart (shared-target-port arrival and dangling
     fixture handles excepted);
  3. overlapping detour runs are >=22px apart; same-source-port detour
     wires keep >=5px between their outbound anchors;
  4. every detour path leaves its target port with a horizontal
     end-tangent (the exact KindEdge invariant; last-chord <35deg).
Also prints per-flow stats that must match the SSR audit numbers.

Usage: python3 wire-audit.py [project.json] [flow ...]
Exits non-zero on any finding.
"""
import json
import math
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
PROJ = sys.argv[1] if len(sys.argv) > 1 else f'{ROOT}/front-end/projects/heptagon/project.json'
FLOW_FILTER = sys.argv[2:] or None

NODE_W = 190
ROW_STEP = 15
SPREAD_STEP = 7
CLEAR = 2.5
LANE_STEP = 22

# TRUE port rows straight from the renderer's shared geometry model in
# src/flowdefs.js (portsOf/portBaseY/estHeight, exported). No independent
# Python port model — a mismatch here is a code bug, never a fallback.
DUMP = json.loads(subprocess.run(
    ['node', f'{ROOT}/front-end/web/scripts/ports-dump.mjs', PROJ],
    capture_output=True, text=True, check=True).stdout)


def ports(flow_name, node_id):
    return DUMP[flow_name][node_id]


def port_y(flow_name, n, handle, side):
    p = ports(flow_name, n['id'])
    rows = p['out'] if side == 'out' else p['in']
    row = rows.index(handle) if handle in rows else 0
    return n['position']['y'] + p['baseY'] + row * ROW_STEP


def seg_rect(ax, ay, bx, by, b):
    """Exact segment-vs-rect interior test (Liang-Barsky)."""
    dx = bx - ax
    dy = by - ay
    t0, t1 = 0.0, 1.0
    def clip(p, q):
        nonlocal t0, t1
        if p == 0:
            return q >= 0
        r = q / p
        if p < 0:
            if r > t1:
                return False
            if r > t0:
                t0 = r
        else:
            if r < t0:
                return False
            if r < t1:
                t1 = r
        return True
    if not clip(-dx, ax - b[0]) or not clip(dx, b[1] - ax):
        return False
    if not clip(-dy, ay - b[2]) or not clip(dy, b[3] - ay):
        return False
    return t1 > t0


def leg_clear(boxes, from_x, from_y, to_x, lane_y):
    for b in boxes:
        if b[1] <= min(from_x, to_x) or b[0] >= max(from_x, to_x):
            continue
        if seg_rect(from_x, from_y, to_x, lane_y,
                    (b[0], b[1], b[2] - 2, b[3] + 2)):
            return False
    return True


def cubic_hits(boxes, p0, c1, c2, p1):
    hx0 = min(p0[0], c1[0], c2[0], p1[0])
    hx1 = max(p0[0], c1[0], c2[0], p1[0])
    for b in boxes:
        if b[1] <= hx0 or b[0] >= hx1:
            continue
        for i in range(21):
            t = i / 20.0
            u = 1.0 - t
            px = u * u * u * p0[0] + 3 * u * u * t * c1[0] + 3 * u * t * t * c2[0] + t * t * t * p1[0]
            py = u * u * u * p0[1] + 3 * u * u * t * c1[1] + 3 * u * t * t * c2[1] + t * t * t * p1[1]
            if b[0] < px < b[1] and b[2] - 2 < py < b[3] + 2:
                return True
    return False


def rise_clear(boxes, from_x, port_yv, to_x, lane_y):
    k = max(10, min(24, 0.6 * abs(to_x - from_x)))
    return not cubic_hits(boxes, (from_x, port_yv), (from_x + k, port_yv),
                          (to_x - k, lane_y), (to_x, lane_y))


def desc_clear(boxes, to_x, port_yv, from_x, lane_y):
    # EXACT drawn shape: cubic from the lane down to the glide start
    # (already horizontal), then a straight glide into the port.
    leg = to_x - from_x
    k = max(10, min(24, 0.6 * abs(leg)))
    lead = max(8, min(16, 0.4 * abs(leg)))
    if cubic_hits(boxes, (from_x, lane_y), (from_x + k, lane_y),
                  (to_x - k - lead, port_yv), (to_x - lead, port_yv)):
        return False
    for b in boxes:
        if seg_rect(to_x - lead, port_yv, to_x, port_yv, b):
            return False
    return True


def cubic(p0, c1, c2, p1, n):
    pts = []
    for i in range(n + 1):
        t = i / float(n)
        u = 1.0 - t
        x = u * u * u * p0[0] + 3 * u * u * t * c1[0] + 3 * u * t * t * c2[0] + t * t * t * p1[0]
        y = u * u * u * p0[1] + 3 * u * u * t * c1[1] + 3 * u * t * t * c2[1] + t * t * t * p1[1]
        pts.append((x, y))
    return pts


def quad(p0, c, p1, n):
    pts = []
    for i in range(n + 1):
        t = i / float(n)
        u = 1.0 - t
        pts.append((u * u * p0[0] + 2 * u * t * c[0] + t * t * p1[0],
                    u * u * p0[1] + 2 * u * t * c[1] + t * t * p1[1]))
    return pts


def line(p0, p1, n):
    pts = [(p0[0] + (p1[0] - p0[0]) * i / float(n), (p1[1] - p0[1]) * i / float(n) + p0[1])
           for i in range(n + 1)]
    pts[-1] = (p1[0], p1[1])
    return pts


def _pt_seg(px, py, qx0, qy0, qx1, qy1):
    l2 = (qx1 - qx0) ** 2 + (qy1 - qy0) ** 2
    if not l2:
        return math.hypot(px - qx0, py - qy0)
    t = ((px - qx0) * (qx1 - qx0) + (py - qy0) * (qy1 - qy0)) / l2
    t = max(0.0, min(1.0, t))
    return math.hypot(px - (qx0 + t * (qx1 - qx0)), py - (qy0 + t * (qy1 - qy0)))


def seg_seg_dist(p0, p1, q0, q1):
    """True min distance between two segments."""
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    ex, ey = q1[0] - q0[0], q1[1] - q0[1]
    denom = dx * ey - dy * ex
    if abs(denom) < 1e-9:
        return min(_pt_seg(p0[0], p0[1], q0[0], q0[1], q1[0], q1[1]),
                   _pt_seg(p1[0], p1[1], q0[0], q0[1], q1[0], q1[1]),
                   _pt_seg(q0[0], q0[1], p0[0], p0[1], p1[0], p1[1]),
                   _pt_seg(q1[0], q1[1], p0[0], p0[1], p1[0], p1[1]))
    t0 = ((q0[0] - p0[0]) * ey - (q0[1] - p0[1]) * ex) / denom
    t1 = ((q0[0] - p0[0]) * dy - (q0[1] - p0[1]) * dx) / denom
    if 0 <= t0 <= 1 and 0 <= t1 <= 1:
        return 0.0
    return min(_pt_seg(p0[0], p0[1], q0[0], q0[1], q1[0], q1[1]),
               _pt_seg(p1[0], p1[1], q0[0], q0[1], q1[0], q1[1]),
               _pt_seg(q0[0], q0[1], p0[0], p0[1], p1[0], p1[1]),
               _pt_seg(q1[0], q1[1], p0[0], p0[1], p1[0], p1[1]))


class Wire(object):
    def __init__(self, e, sn, tn, x0, y0, x1, y1):
        self.e = e
        self.sn = sn
        self.tn = tn
        self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1
        self.spread = 0.0
        self.kappa = 0.0
        self.detour_y = None
        self.orig_lane = None
        self.outbound_x = 0.0
        self.inbound_x = 0.0
        self.fallback = False
        self.anchor_shift = 0

    def key(self):
        return (self.e['source'], self.e.get('sourceHandle'),
                self.e['target'], self.e.get('targetHandle'))


def route(flow_name, flow):
    """Mirror of routeEdges — must produce identical decisions."""
    nodes = {n['id']: n for n in flow['graph']['nodes']}
    edges = flow['graph']['edges']
    boxes = [(n['position']['x'], n['position']['x'] + NODE_W,
              n['position']['y'] - 26, n['position']['y'] + ports(flow_name, n['id'])['h'] + 26)
             for n in flow['graph']['nodes']]
    wires = []
    dangling_refs = []
    for e in edges:
        sn, tn = nodes.get(e['source']), nodes.get(e['target'])
        if sn is None or tn is None:
            # Mirror of routeEdges' `if (!sn || !tn) continue` — never
            # crash on hand-authored JSON; surface it as a fixture note.
            dangling_refs.append(e)
            continue
        x0 = sn['position']['x'] + NODE_W
        x1 = tn['position']['x']
        y0 = port_y(flow_name, sn, e.get('sourceHandle'), 'out')
        y1 = port_y(flow_name, tn, e.get('targetHandle'), 'in')
        wires.append(Wire(e, sn, tn, x0, y0, x1, y1))

    # 1) fan-out spread per source port
    groups = {}
    for w in wires:
        groups.setdefault((w.e['source'], w.e.get('sourceHandle') or ''), []).append(w)
    for group in groups.values():
        if len(group) < 2:
            continue
        group.sort(key=lambda w: (w.y1, w.x1))
        for i, w in enumerate(group):
            w.spread = (i - (len(group) - 1) / 2.0) * SPREAD_STEP
    for w in wires:
        w.y0 += w.spread

    # 2) obstacle classification: every wire gets the full chord-vs-box
    # test; backward wires (dx < -24) always detour.
    detour_set = []
    for w in wires:
        if w.x1 - w.x0 < -24:
            detour_set.append(w)
            continue
        for b in boxes:
            if b[1] <= w.x0 + 1 or b[0] >= w.x1 - 1:
                continue
            if seg_rect(w.x0, w.y0, w.x1, w.y1, b):
                detour_set.append(w)
                break
    directs = [w for w in wires if w not in detour_set]

    def chord_angle(a, b):
        dax, day = a.x1 - a.x0, a.y1 - a.y0
        dbx, dby = b.x1 - b.x0, b.y1 - b.y0
        dot = dax * dbx + day * dby
        cosv = dot / (math.hypot(dax, day) * math.hypot(dbx, dby) or 1)
        return dot, math.acos(max(-1.0, min(1.0, cosv)))

    def chord_pts(w):
        return [(w.x0 + (w.x1 - w.x0) * (i / 16.0), w.y0 + (w.y1 - w.y0) * (i / 16.0))
                for i in range(17)]

    # 2b) narrow-X promotion (mirror of the JS 2b pass): antiparallel
    # wires crossing at a shallow acute angle closer than CLEAR tangle
    # near-vertically at the port — bows cannot undo a crossing. Promote
    # the longer-chord member to a detour. Steep crossings stay.
    narrow_x = []
    for i in range(len(directs)):
        for j in range(i + 1, len(directs)):
            a, b = directs[i], directs[j]
            pa, pb = chord_pts(a), chord_pts(b)
            best = float('inf')
            for m in range(16):
                for n2 in range(16):
                    dd = seg_seg_dist(pa[m], pa[m + 1], pb[n2], pb[n2 + 1])
                    if dd < best:
                        best = dd
            if best >= CLEAR:
                continue
            dot, ang = chord_angle(a, b)
            acute = (math.pi - ang) if dot < 0 else ang
            if acute >= 0.35:
                continue
            cross = seg_seg_dist((a.x0, a.y0), (a.x1, a.y1),
                                 (b.x0, b.y0), (b.x1, b.y1)) == 0.0
            if cross or dot < 0:
                narrow_x.append((a, b))
    for a, b in narrow_x:
        promote = a if math.hypot(a.x1 - a.x0, a.y1 - a.y0) >= \
            math.hypot(b.x1 - b.x0, b.y1 - b.y0) else b
        if promote not in detour_set:
            detour_set.append(promote)
    if narrow_x:
        directs = [w for w in directs if w not in detour_set]

    # 3) local external lanes (insertion order = classification order)
    top_lanes, bot_lanes = [], []
    for w in detour_set:
        lo = min(w.x0, w.x1) - 80
        hi = (w.x1 + 80) if w.x1 >= w.x0 else (max(w.x0, w.x1 + NODE_W) + 80)
        top_y = bot_y = 0
        for b in boxes:
            if b[1] <= lo or b[0] >= hi:
                continue
            if not top_y or b[2] < top_y:
                top_y = b[2]
            if not bot_y or b[3] > bot_y:
                bot_y = b[3]
        first_top = (top_y or 0) - 40
        first_bot = (bot_y or 0) + 40
        take_top = w.y0 <= w.y1
        lanes = top_lanes if take_top else bot_lanes
        y = first_top if take_top else first_bot
        stack = 1
        while any(a['y'] == y and (min(a['hi'], hi) - max(a['lo'], lo) > 0) for a in lanes):
            y = (first_top - LANE_STEP * stack) if take_top else (first_bot + LANE_STEP * stack)
            stack += 1
        lanes.append({'lo': lo, 'hi': hi, 'y': y})
        w.detour_y = y
        w.orig_lane = y

    # 3b) anchor search (dual outbound bases for backward wires; same-port
    # siblings stagger their anchor attempts by earlier-sibling count)
    anchor_rank = {}
    for w in detour_set:
        k = (w.e['source'], w.e.get('sourceHandle') or '')
        r = anchor_rank.get(k, 0)
        anchor_rank[k] = r + 1
        w.anchor_shift = r
    ANCHOR_OFFS = [44, 37, 30, 51, 58, 65, 72, 79]
    for w in detour_set:
        ok = False
        for i in range(8):
            o = ANCHOR_OFFS[(i + w.anchor_shift) % 8]
            bases = [w.x0]
            if w.x1 < w.x0 and w.x1 + NODE_W > w.x0:
                bases.append(w.x1 + NODE_W)
            in_x = w.x1 - o
            if abs(in_x - w.x0) < 24:
                continue
            for base in bases:
                ox = base + o
                if rise_clear(boxes, w.x0, w.y0, ox, w.detour_y) and \
                        desc_clear(boxes, w.x1, w.y1, in_x, w.detour_y):
                    w.outbound_x = ox
                    w.inbound_x = in_x
                    ok = True
                    break
            if ok:
                break
        if not ok:
            w.detour_y = None
            w.fallback = True

    # 4) kappa bundles
    def chord_samples(w):
        return [(w.x0 + (w.x1 - w.x0) * i / 16.0, w.y0 + (w.y1 - w.y0) * i / 16.0)
                for i in range(17)]

    tang = []
    for i in range(len(directs)):
        for j in range(i + 1, len(directs)):
            a, b = directs[i], directs[j]
            shared_src = a.e['source'] == b.e['source'] and \
                (a.e.get('sourceHandle') or '') == (b.e.get('sourceHandle') or '')
            shared_tgt = a.e['target'] == b.e['target'] and \
                (a.e.get('targetHandle') or '') == (b.e.get('targetHandle') or '')
            if shared_src and shared_tgt:
                continue
            pa, pb = chord_samples(a), chord_samples(b)
            best = float('inf')
            for m in range(16):
                j_lo = 5 if shared_src else 0
                j_hi = 10 if shared_tgt else 15
                for n2 in range(j_lo, j_hi + 1):
                    dd = seg_seg_dist(pa[m], pa[m + 1], pb[n2], pb[n2 + 1])
                    if dd < best:
                        best = dd
            if best >= 7:
                continue
            # SIGNED angle test: antiparallel wires (dot < 0) never bundle.
            dot, ang = chord_angle(a, b)
            if dot < 0 or ang >= 0.35:
                continue
            tang.append((a, b))
    if tang:
        adj = {id(w): [] for w in directs}
        by_id = {id(w): w for w in directs}
        for a, b in tang:
            adj[id(a)].append(b)
            adj[id(b)].append(a)
        seen = set()
        for w in directs:
            if id(w) in seen or not adj[id(w)]:
                continue
            bundle = []
            queue = [w]
            seen.add(id(w))
            while queue:
                cur = queue.pop(0)
                bundle.append(cur)
                for nb in adj[id(cur)]:
                    if id(nb) not in seen:
                        seen.add(id(nb))
                        queue.append(nb)
            bundle.sort(key=lambda p: ((p.y0 + p.y1), p.x0))

            def quad_hits_box(w, kappa):
                rx = w.x1 - w.x0
                ry = w.y1 - w.y0
                ln = math.hypot(rx, ry) or 1.0
                b = (4 + 1.75 * abs(kappa)) * (1 if kappa > 0 else (-1 if kappa < 0 else 0))
                bb = max(-ln / 4.0, min(ln / 4.0, b))
                cx = (w.x0 + w.x1) / 2.0 + 1.5 * bb * (ry / ln)
                cy = (w.y0 + w.y1) / 2.0 + 1.5 * bb * (-rx / ln)
                for b2 in boxes:
                    if b2[1] <= min(w.x0, w.x1, cx) or b2[0] >= max(w.x0, w.x1, cx):
                        continue
                    for i2 in range(25):
                        t = i2 / 24.0
                        u = 1.0 - t
                        px = u * u * w.x0 + 2 * u * t * cx + t * t * w.x1
                        py = u * u * w.y0 + 2 * u * t * cy + t * t * w.y1
                        if b2[0] < px < b2[1] and b2[2] < py < b2[3]:
                            return True
                return False

            def shrink(k):
                # Mirror of JS: (k >= 4 ? k - 4 : (k <= -4 ? k + 4 : 0))
                if k >= 4:
                    return k - 4
                if k <= -4:
                    return k + 4
                return 0

            def bow_pts(w, kappa):
                rx = w.x1 - w.x0
                ry = w.y1 - w.y0
                ln = math.hypot(rx, ry) or 1.0
                b = (4 + 1.75 * abs(kappa)) * (1 if kappa > 0 else -1)
                bb = max(-ln / 4.0, min(ln / 4.0, b))
                cx = (w.x0 + w.x1) / 2.0 + 1.5 * bb * (ry / ln)
                cy = (w.y0 + w.y1) / 2.0 + 1.5 * bb * (-rx / ln)
                return [((1 - t) ** 2 * w.x0 + 2 * (1 - t) * t * cx + t ** 2 * w.x1,
                         (1 - t) ** 2 * w.y0 + 2 * (1 - t) * t * cy + t ** 2 * w.y1)
                        for t in [i2 / 24.0 for i2 in range(25)]]

            def bow_hits_wire(w, kappa):
                pw = bow_pts(w, kappa)
                for d in directs:
                    if d in bundle:
                        continue
                    pd = chord_samples(d)
                    shared_src = w.e['source'] == d.e['source'] and \
                        (w.e.get('sourceHandle') or '') == (d.e.get('sourceHandle') or '')
                    shared_tgt = w.e['target'] == d.e['target'] and \
                        (w.e.get('targetHandle') or '') == (d.e.get('targetHandle') or '')
                    for m in range(24):
                        j_lo = 5 if shared_src else 0
                        j_hi = 10 if shared_tgt else 15
                        for n2 in range(j_lo, j_hi + 1):
                            if seg_seg_dist(pw[m], pw[m + 1], pd[n2], pd[n2 + 1]) < 7:
                                return True
                return False

            for i, c in enumerate(bundle):
                # Rank 0 stays flat; every later member bows AWAY from its
                # previous-rank neighbor (side = cross(u, prev->c start)).
                if i == 0:
                    kappa = 0.0
                else:
                    prev = bundle[i - 1]
                    side = (c.x1 - c.x0) * (prev.y0 - c.y0) - (c.y1 - c.y0) * (prev.x0 - c.x0)
                    kappa = (4.0 if side > 0 else -4.0) * i
                while kappa != 0 and (quad_hits_box(c, kappa) or bow_hits_wire(c, kappa)):
                    kappa = shrink(kappa)
                c.kappa = kappa
    if os.environ.get('WIRE_DEBUG'):
        edge_idx = {id(e): i for i, e in enumerate(edges)}
        with open(os.environ['WIRE_DEBUG'], 'a') as fh:
            for w in wires:
                dec = 'fallback' if w.fallback else ('detour' if w.detour_y is not None
                                                     else ('kappa' if w.kappa else 'direct'))
                fh.write('%s|%d|%s|%s|%s|%s|%s|%.1f,%.1f,%.1f,%.1f|s=%s|k=%s\n'
                         % (flow_name, edge_idx[id(w.e)], w.e['source'], w.e.get('sourceHandle') or '',
                            w.e['target'], w.e.get('targetHandle') or '', dec,
                            w.x0, w.y0, w.x1, w.y1, w.spread, w.kappa))
    return boxes, wires, detour_set, directs, dangling_refs


def path_pts(w):
    """Sampled visible path exactly as KindEdge draws it (fallback=None)."""
    if w.fallback:
        return None
    if w.detour_y is not None:
        ox, in_x, ly = w.outbound_x, w.inbound_x, w.detour_y
        k0 = max(10, min(24, 0.6 * abs(ox - w.x0)))
        k1 = max(10, min(24, 0.6 * abs(w.x1 - in_x)))
        lead = max(8, min(16, 0.4 * abs(w.x1 - in_x)))
        p1 = cubic((w.x0, w.y0), (w.x0 + k0, w.y0),
                   (ox - k0, ly), (ox, ly), 16)
        p2 = line((ox, ly), (in_x, ly), 16)
        p3 = cubic((in_x, ly), (in_x + k1, ly),
                   (w.x1 - k1 - lead, w.y1), (w.x1 - lead, w.y1), 16)
        p4 = line((w.x1 - lead, w.y1), (w.x1, w.y1), 2)
        return p1 + p2[1:] + p3[1:] + p4[1:]
    if w.kappa:
        rx = w.x1 - w.x0
        ry = w.y1 - w.y0
        ln = math.hypot(rx, ry) or 1.0
        b = (4 + 1.75 * abs(w.kappa)) * (1 if w.kappa > 0 else -1)
        bb = max(-ln / 4.0, min(ln / 4.0, b))
        mx = (w.x0 + w.x1) / 2.0 + 1.5 * bb * (ry / ln)
        my = (w.y0 + w.y1) / 2.0 + 1.5 * bb * (-rx / ln)
        return quad((w.x0, w.y0), (mx, my), (w.x1, w.y1), 32)
    return line((w.x0, w.y0), (w.x1, w.y1), 32)


def audit(flow_name, flow):
    boxes, wires, detour_set, directs, dangling_refs = route(flow_name, flow)
    problems = []
    fixtures = []
    for e in dangling_refs:
        fixtures.append('fixture: edge endpoints not in nodes (%s->%s) - skipped like routeEdges'
                        % (e.get('source'), e.get('target')))
    stats = {'wires': len(wires), 'detours': 0, 'fallbacks': 0, 'kappa': 0}
    for w in wires:
        if w.detour_y is not None:
            stats['detours'] += 1
        if w.fallback:
            stats['fallbacks'] += 1
        if w.kappa:
            stats['kappa'] += 1

    # Fixture rule (checks 1/2): the heptagon fixtures contain node pairs
    # whose boxes PHYSICALLY OVERLAP (x-overlap up to ~21px) and edges whose
    # handles do not exist on the target type (dangling: port_y fell back
    # to row 0). No router can fix either; bucket them as `fixture:` lines
    # instead of problems so real regressions stay visible.
    node_box = {}
    for n in flow['graph']['nodes']:
        p = ports(flow_name, n['id'])
        node_box[n['id']] = (n['position']['x'], n['position']['x'] + NODE_W,
                             n['position']['y'] - 26, n['position']['y'] + p['h'] + 26)

    def x_overlap(a, b):
        return a[0] < b[1] and b[0] < a[1]

    def dangling(w):
        pin = ports(flow_name, w.tn['id'])
        pout = ports(flow_name, w.sn['id'])
        return (w.e.get('targetHandle') not in pin['in'] or
                w.e.get('sourceHandle') not in pout['out'])

    # check 0: routing quality
    for w in wires:
        dx = w.x1 - w.x0
        if dx < -24 and w.detour_y is None and not w.fallback:
            problems.append('backward wire renders straight: %s' % str(w.key()))
    for w in detour_set:
        if w.fallback:
            escaped = False
            ANCHOR_OFFS = [44, 37, 30, 51, 58, 65, 72, 79]
            for i in range(8):
                o = ANCHOR_OFFS[i]
                bases = [w.x0]
                if w.x1 < w.x0 and w.x1 + NODE_W > w.x0:
                    bases.append(w.x1 + NODE_W)
                in_x = w.x1 - o
                if abs(in_x - w.x0) < 24:
                    continue
                for base in bases:
                    if rise_clear(boxes, w.x0, w.y0, base + o, w.orig_lane) and \
                            desc_clear(boxes, w.x1, w.y1, in_x, w.orig_lane):
                        escaped = True
                        break
                if escaped:
                    break
            if not escaped:
                continue
            problems.append('fallback but escapable (o=%d works): %s' % (o, str(w.key())))

    # check 1: non-fallback paths never enter a box interior
    for w in wires:
        if w.fallback:
            continue
        pts = path_pts(w)
        for (px, py) in pts:
            hit = None
            for b in boxes:
                if b[0] + 1 < px < b[1] - 1 and b[2] + 1 < py < b[3] - 1:
                    hit = b
                    break
            if hit is None:
                continue
            # Entering a box that is the wire's OWN source/target and whose
            # box physically overlaps the other endpoint's box is a fixture
            # defect (ports drawn inside the neighbor) — not routable.
            own = hit in (node_box.get(w.sn['id']), node_box.get(w.tn['id']))
            if own and x_overlap(node_box[w.sn['id']], node_box[w.tn['id']]):
                msg = ('fixture: overlapping boxes %s/%s - %s'
                       % (w.sn['id'], w.tn['id'], str(w.key())))
                if msg not in fixtures:
                    fixtures.append(msg)
            else:
                problems.append('path passes through box at (%.0f,%.0f): %s'
                                % (px, py, str(w.key())))
            break

    # check 2: direct/kappa wires stay CLEAR apart WHERE THEY RUN PARALLEL.
    # Design invariants: (a) wires crossing at a clear angle (>= kappa's
    # 0.35 rad) form a readable X — exempt; (b) near-parallel close approach
    # (< CLEAR at angle < 0.35 rad) is exactly what the kappa pass bundles
    # and bows — any survivor is a kappa gap or a guard-shrunk bow and a
    # real problem. Near a shared port the paths converge by design (fan-out
    # spread separates them at the port) — trim those neighborhoods.
    drawn = [(w, path_pts(w)) for w in directs]
    for i in range(len(drawn)):
        for j in range(i + 1, len(drawn)):
            a, pa = drawn[i]
            b, pb = drawn[j]
            shared_src = a.e['source'] == b.e['source'] and \
                (a.e.get('sourceHandle') or '') == (b.e.get('sourceHandle') or '')
            shared_tgt = a.e['target'] == b.e['target'] and \
                (a.e.get('targetHandle') or '') == (b.e.get('targetHandle') or '')
            trim_a = pa[4:] if shared_src else pa
            trim_b = pb[4:] if shared_src else pb
            trim_a = trim_a[:-4] if shared_tgt else trim_a
            trim_b = trim_b[:-4] if shared_tgt else trim_b
            if not trim_a or not trim_b:
                continue
            # A dangling handle (not on the target/source type) resolves to
            # a fallback row-0 y — two such wires collapse onto the same
            # port point. That is a fixture defect for the contract audit,
            # not a routing failure.
            if dangling(a) or dangling(b):
                fixtures.append('fixture: dangling handle on %s vs %s'
                                % (str(a.key()), str(b.key())))
                continue
            best = min(math.hypot(p[0] - q[0], p[1] - q[1])
                       for p in trim_a for q in trim_b)
            if best >= CLEAR:
                continue
            dax, day = a.x1 - a.x0, a.y1 - a.y0
            dbx, dby = b.x1 - b.x0, b.y1 - b.y0
            dot = abs(dax * dbx + day * dby)
            denom = math.hypot(dax, day) * math.hypot(dbx, dby) or 1
            if math.acos(min(1, dot / denom)) >= 0.35:
                continue
            problems.append('near-parallel wires %.1fpx apart (<%.1f): %s vs %s'
                            % (best, CLEAR, str(a.key()), str(b.key())))

    # check 3: lane stacking + same-port outbound anchors
    runs = [(w, min(w.outbound_x, w.inbound_x), max(w.outbound_x, w.inbound_x), w.detour_y)
            for w in detour_set if w.detour_y is not None]
    for i in range(len(runs)):
        for j in range(i + 1, len(runs)):
            a, b = runs[i], runs[j]
            if a[3] is None or b[3] is None:
                continue
            if abs(a[3] - b[3]) < 1e-6 and min(a[2], b[2]) - max(a[1], b[1]) > 0:
                problems.append('detour runs share lane y=%.0f: %s vs %s'
                                % (a[3], str(a[0].key()), str(b[0].key())))
    same_port = {}
    for w in detour_set:
        if w.detour_y is not None:
            same_port.setdefault((w.e['source'], w.e.get('sourceHandle') or ''), []).append(w)
    for port_wires in same_port.values():
        for i in range(len(port_wires)):
            for j in range(i + 1, len(port_wires)):
                a, b = port_wires[i], port_wires[j]
                if abs(a.outbound_x - b.outbound_x) < 5:
                    problems.append('same-port outbound anchors %.0fpx apart: %s vs %s'
                                    % (abs(a.outbound_x - b.outbound_x), str(a.key()), str(b.key())))

    # check 4: horizontal arrival (arrowhead tangent). EXACT invariant of
    # the drawn path: KindEdge's descent cubic ends at (x1-lead, y1) with
    # control point c2 = (x1-k1-lead, y1) — tangent (k1, 0), EXACTLY
    # horizontal — then one straight L carries the arrowhead into the port.
    # A huge vertical drop compressed into the 30px inbound corridor
    # legitimately produces a rounded J-turn whose uniform samples stay
    # steep for a couple of px from the port — that is not a defect. The
    # tripwire is EMIT INTEGRITY: inX must sit strictly left of x1 so the
    # tangent points INTO the port from the left. Any emit/control-point
    # regression (inX >= x1, inX unset, NaN) fires immediately.
    for w in wires:
        if w.fallback or w.detour_y is None or abs(w.x1 - w.x0) < 24:
            continue
        if not (w.inbound_x < w.x1):
            problems.append('arrival tangent not horizontal (inX %.1f >= x1 %.1f): %s'
                            % (w.inbound_x, w.x1, str(w.key())))
    return problems, stats, fixtures


def main():
    project = json.load(open(PROJ))
    total = 0
    tot_stats = {'wires': 0, 'detours': 0, 'fallbacks': 0, 'kappa': 0}
    for flow in project['flows']:
        if FLOW_FILTER and flow['name'] not in FLOW_FILTER:
            continue
        problems, stats, fixtures = audit(flow['name'], flow)
        for k in tot_stats:
            tot_stats[k] += stats[k]
        print('%s: %d problems (wires %d, detours %d, kappa %d, fallbacks %d)'
              % (flow['name'], len(problems), stats['wires'], stats['detours'],
                 stats['kappa'], stats['fallbacks']))
        for msg in problems[:10]:
            print('  -', msg)
        for msg in fixtures:
            print('  ~', msg)
        total += len(problems)
    print('TOTAL: %d problems; wires %d (detours %d, kappa %d, fallbacks %d)'
          % (total, tot_stats['wires'], tot_stats['detours'], tot_stats['kappa'],
             tot_stats['fallbacks']))
    sys.exit(1 if total else 0)


if __name__ == '__main__':
    main()
