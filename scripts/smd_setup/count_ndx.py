import sys
n = 0
for line in open(sys.argv[1]):
    s = line.strip()
    if s.startswith("[") or not s:
        continue
    n += len(s.split())
print(n)
