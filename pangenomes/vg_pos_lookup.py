import sys
import json
import re

def parse_vg_json(json_file, pos_input_file, ref_path):
    queries = []
    with open(pos_input_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            queries.append({
                'pos':     int(parts[0]),
                'p28_hap': parts[1],
                'p11_hap': parts[2],
            })

    # Parse all JSON objects — may be one large merged subgraph or multiple small ones
    all_paths = []
    with open(json_file) as f:
        content = f.read().strip()
    decoder = json.JSONDecoder()
    idx = 0
    while idx < len(content):
        while idx < len(content) and content[idx] in ' \t\n\r':
            idx += 1
        if idx >= len(content):
            break
        try:
            data, end_idx = decoder.raw_decode(content, idx)
            idx = end_idx
        except json.JSONDecodeError:
            break
        all_paths.extend(data.get('path', []))

    def get_contig(name):
        # Format: SAMPLE#0#CONTIG#subgraph_id[start-end]  or  SAMPLE#0#CONTIG[start-end]
        base = re.sub(r'\[\d+-\d+\]$', '', name)
        parts = base.split('#')
        return parts[2] if len(parts) >= 3 else base

    def mapping_ref_len(mapping):
        """Length of this mapping along the reference (from_length)."""
        return sum(e.get('from_length', 0) for e in mapping.get('edit', []))

    def get_path_offset(name):
        """Extract the #NUMBER offset before [start-end].
        e.g. SAMPLE#0#CONTIG#96318774[...] → 96318774
        REF paths have no such offset (SAMPLE#0#CONTIG[...]) → 0."""
        m = re.search(r'#(\d+)\[', name)
        return int(m.group(1)) if m else 0

    # Build index: node_id -> list of (path_name, node_contig_start_0based)
    # Actual contig coordinate = path_offset + seg_start + cumulative_node_offset
    node_to_paths = {}
    for path in all_paths:
        name = path['name']
        m = re.search(r'\[(\d+)-(\d+)\]$', name)
        if not m:
            continue
        path_offset = get_path_offset(name)
        seg_start   = int(m.group(1))
        cumoff = 0
        for mapping in path.get('mapping', []):
            nid = str(mapping['position']['node_id'])
            node_contig_start = path_offset + seg_start + cumoff
            if nid not in node_to_paths:
                node_to_paths[nid] = []
            node_to_paths[nid].append((name, node_contig_start))
            cumoff += mapping_ref_len(mapping)

    results = []

    for q in queries:
        ref_pos = q['pos']   # 1-based VCF position
        p28_hap = q['p28_hap']
        p11_hap = q['p11_hap']
        target  = ref_pos - 1  # convert to 0-based

        pos28 = None; con28 = None
        pos11 = None; con11 = None

        # Step 1: find the node_id in the REF path that covers `target`
        ref_node_id     = None
        ref_within_node = 0   # 0-based offset of target within the node

        for path in all_paths:
            name = path['name']
            if ref_path not in name:
                continue
            m = re.search(r'\[(\d+)-(\d+)\]$', name)
            if not m:
                continue
            seg_start = int(m.group(1))
            seg_end   = int(m.group(2))
            if not (seg_start <= target < seg_end):
                continue
            # Walk this path segment's mappings to find the exact node
            cumoff = 0
            for mapping in path.get('mapping', []):
                nlen       = mapping_ref_len(mapping)
                node_start = seg_start + cumoff
                if node_start <= target < node_start + nlen:
                    ref_node_id     = str(mapping['position']['node_id'])
                    ref_within_node = target - node_start
                    break
                cumoff += nlen
            break  # stop after first matching REF path segment

        if ref_node_id is None or ref_node_id not in node_to_paths:
            results.append(('.', '.', '.', '.'))
            continue

        # Step 2: find ALT paths that pass through the same node
        for (path_name, node_contig_start) in node_to_paths[ref_node_id]:
            is_p28 = p28_hap in path_name
            is_p11 = p11_hap in path_name
            if not is_p28 and not is_p11:
                continue
            # 1-based position in the ALT contig
            alt_pos = node_contig_start + ref_within_node + 1
            contig  = get_contig(path_name)
            if is_p28 and pos28 is None:
                pos28 = alt_pos
                con28 = contig
            if is_p11 and pos11 is None:
                pos11 = alt_pos
                con11 = contig

        results.append((
            str(pos28) if pos28 is not None else '.',
            con28      if con28 is not None else '.',
            str(pos11) if pos11 is not None else '.',
            con11      if con11 is not None else '.',
        ))

    return results

if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.stderr.write(
            'Usage: python3 vg_pos_lookup.py <vg_json> <pos_input> <ref_path>\n'
        )
        sys.exit(1)

    results = parse_vg_json(sys.argv[1], sys.argv[2], sys.argv[3])
    for pos28, con28, pos11, con11 in results:
        print(f'{pos28}\t{con28}\t{pos11}\t{con11}')
