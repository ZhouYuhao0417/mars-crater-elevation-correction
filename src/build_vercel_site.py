from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FRAGMENT = ROOT / "interactive" / "crater-corrected-3d-paper.html"
INDEX = ROOT / "index.html"
STANDALONE = ROOT / "interactive" / "crater-corrected-3d.html"


fragment = FRAGMENT.read_text(encoding="utf-8")

page = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="Interactive 3D detrended Mars crater terrain model with mapped rim, suspected outflow trace, and projected outlet.">
  <link rel="icon" href="data:,">
  <title>Mars Crater 3D Elevation Correction</title>
  <style>
    :root {{
      color-scheme: dark;
      --background: #0b0d10;
      --foreground: #f3efe7;
      --card: #15191f;
      --card-foreground: #f3efe7;
      --popover: #10141a;
      --popover-foreground: #f3efe7;
      --primary: #d89b57;
      --primary-foreground: #17100a;
      --secondary: #26303a;
      --secondary-foreground: #f3efe7;
      --muted: #34404b;
      --muted-foreground: #aeb7bf;
      --accent: #213342;
      --accent-foreground: #e8f3ff;
      --border: rgba(243,239,231,.18);
      --input: rgba(243,239,231,.22);
      --ring: #d89b57;
      --viz-series-1: #556f88;
      --viz-series-2: #31b7c2;
      --viz-series-3: #d89b57;
      --viz-series-4: #f3d37b;
      --viz-series-5: #f06f5a;
      --viz-series-6: #7ecf8e;
      --spacing-xs: .45rem;
      --spacing-sm: .75rem;
      --spacing-md: 1rem;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--background);
      color: var(--foreground);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      min-height: 100vh;
      background:
        radial-gradient(circle at 18% 12%, rgba(216,155,87,.14), transparent 24rem),
        linear-gradient(135deg, #0b0d10 0%, #111820 48%, #0c1014 100%);
      color: var(--foreground);
    }}
    .page {{
      min-height: 100vh;
      display: grid;
      grid-template-rows: auto 1fr auto;
    }}
    .topbar {{
      display: flex;
      align-items: end;
      justify-content: space-between;
      gap: 1rem;
      padding: 1rem clamp(1rem, 3vw, 2.25rem) .65rem;
      border-bottom: 1px solid var(--border);
      background: rgba(11,13,16,.78);
      backdrop-filter: blur(14px);
    }}
    h1 {{
      margin: 0;
      font-size: clamp(1.05rem, 2vw, 1.55rem);
      font-weight: 500;
      letter-spacing: 0;
    }}
    .subtitle {{
      margin: .25rem 0 0;
      color: var(--muted-foreground);
      font-size: .92rem;
    }}
    .meta {{
      color: var(--muted-foreground);
      font-size: .85rem;
      text-align: right;
      white-space: nowrap;
    }}
    .viewer {{
      padding: .75rem clamp(.5rem, 1.8vw, 1.5rem) 1rem;
    }}
    #crater-paper-3d {{
      min-height: calc(100vh - 9.8rem);
    }}
    #crater-paper-3d #cp3-plot {{
      height: calc(100vh - 13rem);
      min-height: 520px;
    }}
    .viz-controls {{
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: .6rem .75rem;
      padding: .2rem 0 .45rem;
    }}
    .btn {{
      display: inline-flex;
      align-items: center;
      justify-content: center;
      gap: .35rem;
      min-height: 2.25rem;
      border: 1px solid var(--border);
      border-radius: 7px;
      padding: .45rem .7rem;
      background: rgba(243,239,231,.06);
      color: var(--foreground);
      font: inherit;
      cursor: pointer;
    }}
    .btn:hover {{ background: rgba(243,239,231,.1); }}
    .btn-primary {{
      border-color: color-mix(in srgb, var(--primary) 70%, transparent);
      background: var(--primary);
      color: var(--primary-foreground);
    }}
    .btn-ghost {{ background: transparent; }}
    .form-label {{
      display: inline-flex;
      align-items: center;
      gap: .5rem;
      color: var(--muted-foreground);
      font-size: .9rem;
      white-space: nowrap;
    }}
    .form-range {{
      width: 8.5rem;
      accent-color: var(--primary);
    }}
    .text-small {{ font-size: .9rem; }}
    .text-muted {{ color: var(--muted-foreground); }}
    .sr-only {{
      position: absolute;
      width: 1px;
      height: 1px;
      padding: 0;
      margin: -1px;
      overflow: hidden;
      clip: rect(0, 0, 0, 0);
      white-space: nowrap;
      border: 0;
    }}
    footer {{
      display: flex;
      justify-content: space-between;
      gap: 1rem;
      padding: .7rem clamp(1rem, 3vw, 2.25rem);
      border-top: 1px solid var(--border);
      color: var(--muted-foreground);
      font-size: .82rem;
    }}
    a {{ color: var(--foreground); }}
    @media (max-width: 720px) {{
      .topbar {{
        align-items: start;
        flex-direction: column;
      }}
      .meta {{ text-align: left; white-space: normal; }}
      #crater-paper-3d #cp3-plot {{
        height: 68vh;
        min-height: 460px;
      }}
      .form-label {{
        width: 100%;
        justify-content: space-between;
      }}
      .form-range {{ width: min(52vw, 12rem); }}
      footer {{ flex-direction: column; }}
    }}
  </style>
</head>
<body>
  <main class="page">
    <header class="topbar">
      <div>
        <h1>Mars Crater 3D Elevation Correction</h1>
        <p class="subtitle">Irregular rim model, suspected outflow trace, and detrended MOLA topography.</p>
      </div>
      <div class="meta">Projected outlet: 18.106545°N, 76.040457°E</div>
    </header>
    <section class="viewer" aria-label="Interactive 3D crater model">
{fragment}
    </section>
    <footer>
      <span>Data basis: MOLA MEGDR regional topography and user-mapped JMARS vectors.</span>
      <span>Interpretation: compatible with a low-rim outlet; not diagnostic proof of dam breach.</span>
    </footer>
  </main>
</body>
</html>
"""

INDEX.write_text(page, encoding="utf-8")
STANDALONE.write_text(page, encoding="utf-8")
print(INDEX, INDEX.stat().st_size)
print(STANDALONE, STANDALONE.stat().st_size)
