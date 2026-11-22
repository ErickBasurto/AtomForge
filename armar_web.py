from pathlib import Path

carpeta = Path(__file__).resolve().parent
html = (carpeta / "plantilla_atomforge.html").read_text(encoding="utf-8")
json_str = (carpeta / "datos_atomicos.json").read_text(encoding="utf-8")

if "__DATA_JSON__" not in html:
    raise ValueError("No se encontró el marcador __DATA_JSON__")

(carpeta / "index.html").write_text(html.replace("__DATA_JSON__", json_str, 1), encoding="utf-8")
print("Listo: 'index.html' generado correctamente.")