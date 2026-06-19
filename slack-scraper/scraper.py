import os
import re
import time
from pathlib import Path
from datetime import datetime, timezone

from dotenv import load_dotenv
from playwright.sync_api import sync_playwright

load_dotenv()

WORKSPACE_URL = os.environ["SLACK_WORKSPACE_URL"]
EMAIL = os.environ["SLACK_EMAIL"]
PASSWORD = os.environ["SLACK_PASSWORD"]
CHANNEL_URLS = [u.strip() for u in os.environ["CHANNEL_URLS"].split(",") if u.strip()]
HEADLESS = os.environ.get("HEADLESS", "false").lower() == "true"
SCROLL_PASSES = int(os.environ.get("SCROLL_PASSES", "60"))
OUTPUT_DIR = Path("slack-export")

MESSAGE_PANE_SELECTOR = '[data-qa="slack_kit_list"][role="presentation"]'
CHANNEL_NAME_SELECTOR = '[data-qa="channel_name"]'

MESSAGE_ENDPOINTS = ("conversations.history", "conversations.replies")
USER_ENDPOINTS = ("users.info", "users.list", "client.userBoot", "client.boot")


class ChannelCapture:
    """Accumulates messages and user names seen in network traffic for one channel."""

    def __init__(self):
        self.messages = {}
        self.users = {}

    def handle_response(self, response):
        url = response.url
        if not any(ep in url for ep in MESSAGE_ENDPOINTS + USER_ENDPOINTS):
            return
        try:
            data = response.json()
        except Exception:
            return  

        if any(ep in url for ep in MESSAGE_ENDPOINTS):
            msg_list = data.get("messages", [])
            if isinstance(msg_list, dict):
                msg_list = msg_list.get("items", [])
            for msg in msg_list:
                if isinstance(msg, dict) and msg.get("ts"):
                    self.messages[msg["ts"]] = msg 

        if any(ep in url for ep in USER_ENDPOINTS):
            for member in data.get("members", []) or []:
                self._store_user(member)
            user = data.get("user")
            if user:
                self._store_user(user)

    def _store_user(self, member):
        uid = member.get("id")
        profile = member.get("profile", {}) or {}
        name = profile.get("display_name") or profile.get("real_name") or member.get("name")
        if uid and name:
            self.users[uid] = name


def clean_text(text, users):
    """Replace raw <@U123> mentions with resolved display names where possible."""
    if not text:
        return ""

    def repl(match):
        uid = match.group(1)
        return f"@{users.get(uid, uid)}"

    return re.sub(r"<@([A-Z0-9]+)>", repl, text)


def ts_to_iso(ts):
    try:
        return datetime.fromtimestamp(float(ts), tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    except (TypeError, ValueError):
        return ts


def login(page):
    page.goto(f"https://{WORKSPACE_URL}")
    page.wait_for_selector('input[type="email"]', timeout=30000)
    page.fill('input[type="email"]', EMAIL)

    if page.locator('input[type="password"]').count() == 0:
        page.click('button[type="submit"]')
        page.wait_for_selector('input[type="password"]', timeout=30000)

    page.fill('input[type="password"]', PASSWORD)
    page.click('button[type="submit"]')

    print("If Slack shows a CAPTCHA or 2FA prompt, solve it in the browser window now.")
    page.wait_for_selector(MESSAGE_PANE_SELECTOR, timeout=120000)


def scroll_to_top(page):
    pane = page.locator(MESSAGE_PANE_SELECTOR)
    last_height = None
    for _ in range(SCROLL_PASSES):
        pane.evaluate("el => el.scrollTop = 0")
        time.sleep(1.2)  # let Slack lazy-load + fetch older messages
        height = pane.evaluate("el => el.scrollHeight")
        if height == last_height:
            break
        last_height = height


def extract_channel_name(page, fallback):
    try:
        return page.locator(CHANNEL_NAME_SELECTOR).first.inner_text().strip()
    except Exception:
        return fallback


def export_channel(page, url, index):
    capture = ChannelCapture()
    handler = capture.handle_response
    page.on("response", handler)

    try:
        page.goto(url)
        page.wait_for_selector(MESSAGE_PANE_SELECTOR, timeout=30000)
        scroll_to_top(page)
        time.sleep(1.5)  
    finally:
        page.remove_listener("response", handler)

    name = extract_channel_name(page, f"channel-{index}")
    safe_name = "".join(c if c.isalnum() or c in "-_" else "-" for c in name)

    lines = [f"# {name}\n"]
    for ts, msg in sorted(capture.messages.items(), key=lambda kv: float(kv[0])):
        uid = msg.get("user") or msg.get("bot_id")
        user = capture.users.get(msg.get("user"), msg.get("username") or uid or "unknown")
        text = clean_text(msg.get("text"), capture.users)
        is_reply = bool(msg.get("thread_ts") and msg.get("thread_ts") != ts)
        prefix = "    " if is_reply else ""
        lines.append(f"{prefix}**[{ts_to_iso(ts)}] {user}:** {text}")

    OUTPUT_DIR.mkdir(exist_ok=True)
    out_path = OUTPUT_DIR / f"{safe_name}.md"
    out_path.write_text("\n\n".join(lines) + "\n", encoding="utf-8")
    print(f"Saved {out_path} ({len(capture.messages)} messages)")


def main():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=HEADLESS)
        context = browser.new_context()
        page = context.new_page()

        login(page)

        for i, url in enumerate(CHANNEL_URLS, start=1):
            export_channel(page, url, i)

        browser.close()


if __name__ == "__main__":
    main()