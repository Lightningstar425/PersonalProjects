import argparse
import re
from pathlib import Path


SENTENCE_END_RE = re.compile(r'[.!?]["\')\]]*$')


def ends_sentence(word):
    return bool(SENTENCE_END_RE.search(word))


def chunk_words(text, chunk_size, overlap):
    """Build overlapping chunks that are *at least* chunk_size words AND
    always end on real sentence-ending punctuation. Never truncates mid-
    sentence: if the chunk_size-th word doesn't end a sentence, the chunk
    keeps growing word-by-word until one does (or the text runs out)."""
    words = text.split()
    n = len(words)
    if n == 0:
        return []
    if overlap >= chunk_size:
        raise ValueError("overlap must be smaller than chunk_size")

    chunks = []
    start = 0
    while start < n:
        end = min(start + chunk_size, n)
        while end < n and not ends_sentence(words[end - 1]):
            end += 1
        chunks.append(" ".join(words[start:end]))

        if end >= n:
            break

        # measure overlap from the chunk's *actual* end, since it may have
        # grown past chunk_size to finish its sentence
        next_start = end - overlap
        start = next_start if next_start > start else start + 1
    return chunks


def load_body(path):
    """Strip the leading '# Channel Name' heading, if present, and return
    (channel_name, body_text) so the heading doesn't get counted as content."""
    raw = path.read_text(encoding="utf-8")
    match = re.match(r'^#\s*(.+?)\s*\n+', raw)
    if match:
        return match.group(1).strip(), raw[match.end():].strip()
    return path.stem, raw.strip()


def main():
    parser = argparse.ArgumentParser(description="Chunk cleaned Slack export files for RAG/embedding use.")
    parser.add_argument("--input", default="slack-export-clean", help="Folder of cleaned .md files")
    parser.add_argument("--output", default="slack_chunks", help="Folder to write numbered chunk files")
    parser.add_argument("--chunk-size", type=int, default=600, help="Target words per chunk")
    parser.add_argument("--overlap", type=int, default=90, help="Words of overlap between consecutive chunks")
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)

    md_files = sorted(input_dir.glob("*.md"))
    if not md_files:
        print(f"No .md files found in {input_dir}")
        return

    # First pass: build every chunk in memory so the total count is known
    # up front, which lets filenames be zero-padded consistently.
    all_chunks = []  # (channel_name, chunk_text, chunk_index_within_channel)
    for path in md_files:
        channel_name, body = load_body(path)
        pieces = chunk_words(body, args.chunk_size, args.overlap)
        for i, piece in enumerate(pieces, start=1):
            all_chunks.append((channel_name, piece, i))

    if not all_chunks:
        print("No content to chunk.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    width = max(3, len(str(len(all_chunks))))

    for n, (channel_name, piece, idx) in enumerate(all_chunks, start=1):
        out_path = output_dir / f"slack_file_{n:0{width}d}.md"
        out_path.write_text(f"# {channel_name} (part {idx})\n\n{piece}\n", encoding="utf-8")

    print(f"Wrote {len(all_chunks)} chunk file(s) to {output_dir}/ "
          f"(chunk size {args.chunk_size} words, overlap {args.overlap} words)")


if __name__ == "__main__":
    main()