"""Redact local file paths and user names from the Level 3 and Level 4 session transcripts before a public release.
usage: python3 tools/redact_transcripts.py [--write] [--user=LOGINNAME]     (default: dry run, prints what would change)
Replaces Windows and macOS home-directory paths, Cursor project-folder names derived from them, and the operator's login name with
neutral placeholders.  Edit PATTERNS to add anything else.  The original is kept as <name>.orig unless --write is omitted."""
import os, re, sys
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = ["level3/cursor_instructions_execution_process.md", "level4/cursor_instructions_execution_findings.md"]
PATTERNS = [
    (re.compile(r"[A-Za-z]:\\Users\\[^\\\s\"'`)]+(?:\\[^\s\"'`)]*)?"), "<LOCAL_PATH>"),                 # C:\Users\name\...
    (re.compile(r"/(?:Users|home)/[^/\s\"'`)]+(?:/[^\s\"'`)]*)?"), "<LOCAL_PATH>"),                      # /Users/name/...
    (re.compile(r"c-Users-[A-Za-z0-9_-]+"), "<CURSOR_PROJECT>"),                                           # Cursor project folder slug
]
write = "--write" in sys.argv
for a in sys.argv[1:]:                       # --user=NAME adds a login name to redact wherever it still appears
    if a.startswith("--user="): PATTERNS.append((re.compile(re.escape(a[7:]), re.I), "<USER>"))
for rel in FILES:
    p = os.path.join(root, rel); s = open(p, encoding="utf-8").read(); t = s; n = 0
    for pat, rep in PATTERNS:
        t, k = pat.subn(rep, t); n += k
    print(f"{rel}: {n} replacements" + ("" if n else " (nothing to redact)"))
    if write and t != s:
        open(p + ".orig", "w", encoding="utf-8").write(s); open(p, "w", encoding="utf-8").write(t)
        left = [m.group(0) for m in re.finditer(r"Users[\\/]|home/", t)]
        print(f"   written; {len(left)} identifiers remain" + (f": {left[:3]}" if left else ""))
if not write: print("dry run: nothing written (use --write)")
