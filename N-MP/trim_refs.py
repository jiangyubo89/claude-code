"""
Trim uncited references and renumber everything consistently.
"""
import re, sys
sys.stdout.reconfigure(encoding='utf-8')

with open('N-MP/manuscript/Chapter8_Identification_and_Tracking_Tools.md', 'r', encoding='utf-8') as f:
    content = f.read()

# Split content into body and references
ref_marker = '## References'
ref_idx = content.find(ref_marker)
body = content[:ref_idx]
ref_section = content[ref_idx:]

# Extract all cited reference numbers from body
cited_nums = sorted(set(int(x) for x in re.findall(r'\[(\d+)\]', body)))
print(f'Cited references: {len(cited_nums)}')

# Parse reference entries
ref_entries = {}
for m in re.finditer(r'^(\d+)\.\s+(.*?)(?=\n\d+\.\s|\n---|\n##|\Z)', ref_section, re.MULTILINE | re.DOTALL):
    num = int(m.group(1))
    text = m.group(2).strip()
    ref_entries[num] = text

print(f'Total reference entries: {len(ref_entries)}')

# Keep only cited references
kept = {n: ref_entries[n] for n in cited_nums if n in ref_entries}
print(f'Keeping: {len(kept)} references')

# Create old->new number mapping
old_to_new = {}
for new_num, old_num in enumerate(sorted(kept.keys()), 1):
    old_to_new[old_num] = new_num

print(f'Renumbering: {min(old_to_new.keys())}-{max(old_to_new.keys())} -> 1-{len(old_to_new)}')

# Replace all [N] citations in body with new numbers
# Must handle [N,M] and [N-M] patterns too
def replace_cite(match):
    inner = match.group(1)
    # Handle comma-separated: [42,44]
    parts = inner.split(',')
    new_parts = []
    for p in parts:
        p = p.strip()
        if p.isdigit():
            old = int(p)
            if old in old_to_new:
                new_parts.append(str(old_to_new[old]))
            else:
                new_parts.append(p)  # keep as-is if not found
        else:
            new_parts.append(p)
    return '[' + ','.join(new_parts) + ']'

new_body = re.sub(r'\[([\d,\s]+)\]', replace_cite, body)

# Build new reference list
new_ref_lines = [ref_marker + '\n']
for old_num in sorted(kept.keys()):
    new_num = old_to_new[old_num]
    new_ref_lines.append(f'\n{new_num}. {kept[old_num]}\n')

# Reconstruct: find everything after references (Figure Legends etc)
# Find Figure Legends section in ref_section
fig_legend_marker = '## Figure Legends'
fl_idx = ref_section.find(fig_legend_marker)
if fl_idx > 0:
    after_refs = ref_section[fl_idx:]
    # Also renumber citations in figure legends
    after_refs = re.sub(r'\[([\d,\s]+)\]', replace_cite, after_refs)
else:
    after_refs = ''

new_content = new_body + '\n' + ''.join(new_ref_lines) + '\n---\n\n' + after_refs

with open('N-MP/manuscript/Chapter8_Identification_and_Tracking_Tools.md', 'w', encoding='utf-8') as f:
    f.write(new_content)

# Verify
with open('N-MP/manuscript/Chapter8_Identification_and_Tracking_Tools.md', 'r', encoding='utf-8') as f:
    verify = f.read()

new_cited = sorted(set(int(x) for x in re.findall(r'\[(\d+)\]', verify[:verify.find('## References')])))
new_entries = re.findall(r'^(\d+)\.\s', verify[verify.find('## References'):], re.MULTILINE)
print(f'\n=== RESULT ===')
print(f'References in text: {len(new_cited)}, range: {min(new_cited)}-{max(new_cited)}')
print(f'Reference entries: {len(new_entries)}, range: 1-{max(int(x) for x in new_entries)}')

# Check for any cited but missing
entry_set = set(int(x) for x in new_entries)
missing = set(new_cited) - entry_set
if missing:
    print(f'WARNING - cited but missing: {sorted(missing)}')
else:
    print('All cited references have entries. OK.')
