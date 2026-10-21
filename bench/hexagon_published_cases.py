from pathlib import Path
import re,json,random,argparse
from fractions import Fraction
p=argparse.ArgumentParser(description="Extract independent u11 symbol-word checks from Dixon-Liu EZMHVg_u11_lin.txt.")
p.add_argument('--input',required=True,type=Path)
p.add_argument('--output',required=True,type=Path)
a=p.parse_args();r=a.output.resolve();r.mkdir(parents=True,exist_ok=True)
s=a.input.read_text().replace('\\','')
cases=[]
for loop in [1,2,8]:
 start=re.search(r'EZMHVg_u11_u\['+str(loop)+r'\]\s*:=',s).end();end=s.index(':',start)
 body=re.sub(r'\s+','',s[start:end])
 values={}
 for m in re.finditer(r'HPL\(\[([01](?:,[01])*)\],u\)',body):
  bits=m.group(1).split(',')
  if len(bits)!=2*loop:continue
  prefix=body[max(0,m.start()-100):m.start()]
  q=re.search(r'([+-]?\d+(?:/\d+)?)\*$',prefix)
  if q:coefficient=Fraction(q.group(1))
  else:
   sign=prefix[-1:] if prefix else '+'
   if sign not in '+-':raise RuntimeError('Unsupported top-weight HPL coefficient: '+prefix[-40:])
   coefficient=Fraction(-1 if sign=='-' else 1)
  word=''.join(reversed(bits))
  values[word]=values.get(word,Fraction(0))+coefficient*((-1)**bits.count('1'))
 # Include pure logs, both ends of the alphabet and deterministic mixed words.
 words={'0'*(2*loop),'0'*(2*loop-1)+'1','1'*(2*loop),'01'*loop,'0'+'1'*(2*loop-1)}
 rng=random.Random(20261006+loop)
 keys=sorted(values)
 for word in rng.sample(keys,min(20,len(keys))):words.add(word)
 for word in sorted(words):cases.append({'loop':loop,'word':word,'expected':str(values.get(word,Fraction(0)))})
 print(loop,'top_weight_HPL_terms',len(values),'validation_words',len(words))
(r/'cases.json').write_text(json.dumps(cases,indent=2)+'\n')
with (r/'cases.tsv').open('w') as out:
 for x in cases:out.write(f"{x['loop']}\t{x['word']}\t{x['expected']}\n")
