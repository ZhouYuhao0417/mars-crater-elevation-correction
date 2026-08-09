from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "index.html"
PLOTLY = ROOT / "vendor" / "plotly.min.js"
OUT = ROOT / "deploy" / "vercel-import.html"


page = INDEX.read_text(encoding="utf-8")
plotly = PLOTLY.read_text(encoding="utf-8").replace("</script>", "<\\/script>")
bundle = page.replace(
    '<script src="/vendor/plotly.min.js"></script>',
    f"<script>\n{plotly}\n</script>",
)

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(bundle, encoding="utf-8")
print(OUT, OUT.stat().st_size)
