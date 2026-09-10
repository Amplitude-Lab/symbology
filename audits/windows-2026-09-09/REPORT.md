# Windows compatibility follow-up — 9 September 2026

macOS compatibility is confirmed by the project owner. This follow-up focuses
on Windows. **Windows acceptance remains pending:** the available host runs
Linux and has no Windows, WSL or Windows emulation runtime. The Windows CI job
added in this pass has not been run remotely.

## Scope and corrections

The documented complete Windows workflow runs the core, server and exported
scripts inside WSL, with the webpage opened in a Windows browser. The native
Windows launcher is an editor/server entry point; it does not provide a
complete native port of the numerical tools.

The source audit found and corrected these issues:

- Project names such as `NUL`, `CON` and `COM1` became reserved Windows device
  names. Eight new project/export cases failed before the correction. New
  project directories and export filenames now receive a safe prefix when
  necessary; the user-visible project name is preserved. Existing project
  directories are not renamed. Windows reserves these names even when an
  extension follows them ([Microsoft filename rules](https://learn.microsoft.com/en-us/windows/win32/fileio/naming-a-file)).
- Exported Bash scripts used platform-default text output, which allows CRLF
  on Windows. Exports now explicitly use UTF-8 and LF; embedded Python source
  is read as UTF-8. This fixes script encoding, not arbitrary cross-OS plan
  relocation ([Python text-file API](https://docs.python.org/3.12/library/pathlib.html#pathlib.Path.write_text)).
- The Windows launcher continued after dependency/build errors and treated
  any existing virtual-environment or web-build directory as complete. It
  now checks the Python executable and built index, retries dependencies,
  propagates failures, uses `npm ci`, preserves the server exit code and
  supports suppressing automatic browser opening during tests. Native
  execution of these batch-file changes still awaits Windows CI.

No numerical C++ files or solver algorithms were changed in this follow-up.

## Verification and prepared Windows coverage

| Check | Result |
|---|---|
| Linux API, export and portability suite | 58 passed; 6 native-Windows launcher cases skipped |
| Browser-test runner, real local API and Chromium | 14/14 passed on Linux |
| Python compilation, workflow YAML parsing and source whitespace | Passed |
| Native Windows launcher/API/web build/Chromium | CI job prepared; not run |
| Full WSL build, numerical suite and Windows-browser access | Not run |

The new `windows-editor` CI job uses Windows Server 2022, Python 3.12 and
Node.js 20. It builds the web bundle and runs the portable API tests, native
launcher cases and existing Chromium editor scenarios. Launcher cases cover
paths containing spaces, parentheses, ampersands and apostrophes; incomplete
setup; dependency/build failures; and server exit status. A separate process
must be excluded by the project lock and admitted after release.

The three existing Bash-dependent API cases are skipped in native Windows.
They remain active on Linux/WSL. The Windows browser checks exercise the
editor, draft recovery and saving, not native C++ calculations. The new
cross-platform browser runner starts an isolated server, chooses a local
port, and cleans up the temporary project store and server.

Evidence: [API results](api-results.txt),
[local browser results](browser-local-results.json),
[local browser console](browser-local-results.txt).
Reproduction: [test instructions](../../tests/README.md).

## Remaining Windows limitations and acceptance steps

1. Run the new Windows CI job after these changes are committed and pushed.
   Its configuration is checked locally, but it is not a substitute for a
   successful result on a Windows runner.
2. In a fresh WSL checkout, install the documented Linux dependencies and
   run `make check-public`, the full Python suite and the browser checks.
   Keep the repository and project files in the Linux home directory and
   create the Python and Node environments inside WSL. This follows
   [Microsoft's filesystem guidance](https://learn.microsoft.com/en-us/windows/wsl/filesystems).
   Verify that a Windows browser can open the WSL server, edit/save a flow,
   run the E6 example, cancel a run and execute a relocated export.
3. Generate exports inside WSL. A dry-run probe confirmed that a plan with
   native Windows `C:\...` paths retains those paths when run on POSIX:
   the current relocator assumes the source and destination path separator
   agree. See the [probe result](native-export-path-probe.txt). Native
   Windows-to-Linux export relocation is not certified by the LF fix.
4. A full native Windows numerical port needs separate work: `compute_rhs`
   uses POSIX shell quoting and looks for an unsuffixed sibling `bootstrap`,
   while native Windows process cancellation currently kills only the
   immediate child. These paths are reasons to retain WSL as the supported
   execution route. Windows `system()` uses the command interpreter
   ([Microsoft runtime documentation](https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/system-wsystem)).

The earlier macOS follow-up item is closed based on the owner's confirmation.
This Windows follow-up builds on robustness commit `da66872`.
