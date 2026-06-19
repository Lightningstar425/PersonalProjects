
import argparse
import re
from pathlib import Path

MESSAGE_BLOCK_RE = re.compile(
    r'^(?P<indent>[ \t]*)\*\*\[(?P<ts>[^\]]+)\]\s+(?P<user>.+?):\*\*[ \t]?'
    r'(?P<text>.*?)(?=\n[ \t]*\*\*\[|\Z)',
    re.MULTILINE | re.DOTALL,
)

LINK_RE = re.compile(r'<(https?://[^>|]+)(?:\|([^>]*))?>')
CHANNEL_REF_RE = re.compile(r'<#[A-Z0-9]+(?:\|([^>]*))?>')
USER_REF_RE = re.compile(r'<@[A-Z0-9]+(?:\|([^>]*))?>')
SPECIAL_MENTION_RE = re.compile(r'<![a-zA-Z]+(?:\^[^|>]+)?(?:\|([^>]*))?>')
EMOJI_RE = re.compile(r':[a-zA-Z0-9_+\-]+:')
CODE_BLOCK_RE = re.compile(r'```(.*?)```', re.DOTALL)
INLINE_CODE_RE = re.compile(r'`([^`\n]+)`')
BOLD_RE = re.compile(r'\*([^*\n]+)\*')
ITALIC_RE = re.compile(r'(?<!\w)_([^_\n]+)_(?!\w)')
STRIKE_RE = re.compile(r'~([^~\n]+)~')


def extract_links(text, ts):
    """Pull (timestamp, label, url) tuples out of text, removing the raw
    link markup but keeping any human-readable label inline so the
    sentence still reads naturally."""
    found = []

    def repl(m):
        url, label = m.group(1), m.group(2)
        found.append((ts, label.strip() if label else None, url.strip()))
        return label.strip() if label else ""

    text = LINK_RE.sub(repl, text)
    return text, found


def strip_symbols(text):
    text = CHANNEL_REF_RE.sub(lambda m: f"#{m.group(1).strip()}" if m.group(1) else "", text)
    text = USER_REF_RE.sub(lambda m: m.group(1).strip() if m.group(1) else "", text)
    text = SPECIAL_MENTION_RE.sub(lambda m: m.group(1).strip() if m.group(1) else "", text)
    text = EMOJI_RE.sub("", text)
    text = CODE_BLOCK_RE.sub(lambda m: m.group(1).strip(), text)
    text = INLINE_CODE_RE.sub(lambda m: m.group(1), text)
    text = BOLD_RE.sub(lambda m: m.group(1), text)
    text = ITALIC_RE.sub(lambda m: m.group(1), text)
    text = STRIKE_RE.sub(lambda m: m.group(1), text)
    text = re.sub(r'[ \t]{2,}', ' ', text)        # collapse extra spaces left by removed tags
    text = re.sub(r'[ \t]+\n', '\n', text)        # trailing spaces on a line
    text = re.sub(r'\n{3,}', '\n\n', text)        # collapse big gaps
    return text.strip()


def process_file(path):
    raw = path.read_text(encoding="utf-8")

    heading = ""
    body = raw
    head_match = re.match(r'^(#[^\n]*\n)', raw)
    if head_match:
        heading = head_match.group(1)
        body = raw[head_match.end():]

    cleaned_blocks = []
    all_links = []

    for m in MESSAGE_BLOCK_RE.finditer(body):
        indent = m.group("indent") or ""
        ts = m.group("ts").strip()
        text = m.group("text")

        text, links = extract_links(text, ts)
        all_links.extend(links)

        text = strip_symbols(text)
        if not text:
            continue

        cleaned_blocks.append(f"{indent}{text}")

    cleaned = (heading + "\n" if heading else "") + "\n\n".join(cleaned_blocks) + "\n"
    return cleaned, all_links


def write_links_file(path, links):
    if not links:
        return
    lines = [f"# Links from {path.stem}\n"]
    for ts, label, url in links:
        if label:
            lines.append(f"- {label}: {url}")
        else:
            lines.append(f"- {url}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Clean Slack export Markdown files and extract links.")
    parser.add_argument("--input", default="slack-export", help="Folder containing the raw exported .md files")
    parser.add_argument("--output", default="slack-export-clean", help="Folder to write cleaned .md files")
    parser.add_argument("--links-dir", default="slack_links", help="Folder to write extracted-link files")
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    links_dir = Path(args.links_dir)

    md_files = sorted(input_dir.glob("*.md"))
    if not md_files:
        print(f"No .md files found in {input_dir}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    links_dir.mkdir(parents=True, exist_ok=True)

    for path in md_files:
        cleaned, links = process_file(path)

        out_path = output_dir / path.name
        out_path.write_text(cleaned, encoding="utf-8")

        links_path = links_dir / f"{path.stem}-links.md"
        write_links_file(links_path, links)

        print(f"{path.name}: cleaned -> {out_path}  |  {len(links)} link(s) -> "
              f"{links_path if links else '(none)'}")


if __name__ == "__main__":
    main()