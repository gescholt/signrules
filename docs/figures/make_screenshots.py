"""
Screenshots of the signrules app for the README.

    python docs/figures/make_screenshots.py      # writes docs/figures/app_A-e.png and app_A-e_verdict.png

The app (``signrules.editor.PatternEditor``) is run headless on the preset A-e; the first picture is the
pattern as loaded (sign matrix, network, status line), the second the verdict after ``r``.  The pictures are
Textual's own SVG screenshots, so nothing is drawn by hand.  They are converted to PNG here (needs
``cairosvg``, and the fonts DejaVu Sans Mono and FreeMono), because an SVG on a web page is drawn with the
reader's fonts: the box characters then drift and the braille dots of the network may be missing.
"""

import asyncio
import re
from pathlib import Path

from signrules.editor import PatternEditor
from signrules.patterns import ELEVEN, short_id

OUT = Path(__file__).resolve().parent
NAME = "A-e competitor + 3-cycle (Thm 3.5)"
MONO, BRAILLE = "DejaVu Sans Mono", "FreeMono"          # DejaVu Sans Mono has no braille patterns
WIDTH = 1600                                            # pixels


def to_png(svg: str, stem: str) -> None:
    import cairosvg
    svg = re.sub(r"@font-face\s*{[^}]*}", "", svg)                       # the web fonts Textual points to
    svg = svg.replace('"Fira Code"', f'"{MONO}"').replace("Fira Code, monospace", f'"{MONO}", monospace')
    svg = re.sub(r"<text ([^>]*)>([^<]*[⠀-⣿][^<]*)</text>",
                 lambda m: f'<text font-family="{BRAILLE}" {m.group(1)}>{m.group(2)}</text>', svg)
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(OUT / f"{stem}.png"), output_width=WIDTH)
    print("wrote", stem)


async def shot(stem: str, size, *keys, wait: float = 0.5) -> None:
    a, B = ELEVEN[NAME]
    app = PatternEditor(list(a), [list(map(int, r)) for r in B], short_id(NAME))
    async with app.run_test(size=size) as pilot:
        await pilot.pause(0.5)
        for k in keys:
            await pilot.press(k)
        await pilot.pause(wait)
        to_png(app.export_screenshot(), stem)


async def main() -> None:
    await shot("app_A-e", (104, 28))
    await shot("app_A-e_verdict", (104, 52), "r", wait=3.0)   # the search of A-e ends at the root at once


if __name__ == "__main__":
    asyncio.run(main())
