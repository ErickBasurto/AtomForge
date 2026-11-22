import json
import math
from pathlib import Path
import numpy as np


def generar_sistemas_fcc():
    # Familia {111}<110>: 4 planos compactos x 3 direcciones compactas = 12 sistemas
    planos = [
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ]
    sistemas = []
    for n in planos:
        # Direcciones <110> ortogonales al plano (n · s = 0)
        candidatos = [
            [1, -1, 0], [1, 0, -1], [0, 1, -1],
            [1, 1, 0], [1, 0, 1], [0, 1, 1],
        ]
        dirs_validas = [s for s in candidatos if np.dot(n, s) == 0]
        for s in dirs_validas:
            sistemas.append({"plano": n, "direccion": s})
    return sistemas


def generar_sistemas_bcc():
    # Familia primaria {110}<111>: 6 planos x 2 direcciones = 12 sistemas
    planos = [
        [1, 1, 0], [1, -1, 0],
        [1, 0, 1], [1, 0, -1],
        [0, 1, 1], [0, 1, -1],
    ]
    candidatos = [
        [1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1],
    ]
    sistemas = []
    for n in planos:
        dirs_validas = [s for s in candidatos if np.dot(n, s) == 0]
        for s in dirs_validas:
            sistemas.append({"plano": n, "direccion": s})
    return sistemas


def calcular_schmid_direccion(uvw, sistemas):
    u = np.array(uvw, dtype=float)
    u_norm = np.linalg.norm(u)
    resultados = []
    for idx, sys in enumerate(sistemas, start=1):
        n = np.array(sys["plano"], dtype=float)
        s = np.array(sys["direccion"], dtype=float)
        cos_phi = abs(float(np.dot(u, n))) / (u_norm * float(np.linalg.norm(n)))
        cos_lam = abs(float(np.dot(u, s))) / (u_norm * float(np.linalg.norm(s)))
        m = cos_phi * cos_lam
        resultados.append({
            "id": idx,
            "plano": sys["plano"],
            "direccion": sys["direccion"],
            "cos_phi": round(cos_phi, 4),
            "cos_lambda": round(cos_lam, 4),
            "m": round(m, 4),
        })
    resultados.sort(key=lambda x: x["m"], reverse=True)
    return resultados


def main():
    carpeta = Path(__file__).resolve().parent
    salida = carpeta / "datos_atomicos.json"

    sistemas_fcc = generar_sistemas_fcc()
    sistemas_bcc = generar_sistemas_bcc()

    # 1. Metales de referencia para la Pestaña 1 y simuladores (fuente: Callister / literatura estándar)
    metales_base = [
        {
            "simbolo": "Cu",
            "nombre": "Cobre puro (monocristal)",
            "estructura": "FCC",
            "parametro_red_nm": 0.3615,
            "fea": 0.74,
            "tau_crss_mpa": 0.70,
            "uso": "Conductores eléctricos y estudio clásico de deslizamiento {111}<110>",
        },
        {
            "simbolo": "Al",
            "nombre": "Aluminio puro (monocristal)",
            "estructura": "FCC",
            "parametro_red_nm": 0.4049,
            "fea": 0.74,
            "tau_crss_mpa": 0.78,
            "uso": "Estructuras aeroespaciales ligeras y ductilidad criogénica",
        },
        {
            "simbolo": "Ni",
            "nombre": "Níquel / Superaleación base Ni",
            "estructura": "FCC",
            "parametro_red_nm": 0.3524,
            "fea": 0.74,
            "tau_crss_mpa": 3.20,
            "uso": "Álabes monocristalinos orientados en [001] para turbinas de avión",
        },
        {
            "simbolo": "Fe-α",
            "nombre": "Hierro Ferrítico (a temperatura ambiente)",
            "estructura": "BCC",
            "parametro_red_nm": 0.2866,
            "fea": 0.68,
            "tau_crss_mpa": 27.50,
            "uso": "Matriz base de aceros al carbono; presenta transición dúctil-frágil",
        },
        {
            "simbolo": "Fe-γ",
            "nombre": "Hierro Austenítico (> 912 °C)",
            "estructura": "FCC",
            "parametro_red_nm": 0.3589,
            "fea": 0.74,
            "tau_crss_mpa": 12.00,
            "uso": "Fase huésped para tratamientos de cementación y temple de aceros",
        },
    ]

    # 2. Direcciones precargadas de carga uniaxial (incluyendo casos de tus apuntes)
    direcciones_preset = [
        {"etiqueta": "[1 0 0] — Eje cúbico (orientación de álabes monocristalinos)", "uvw": [1, 0, 0]},
        {"etiqueta": "[1 1 0] — Diagonal de cara", "uvw": [1, 1, 0]},
        {"etiqueta": "[1 1 1] — Diagonal principal del cubo", "uvw": [1, 1, 1]},
        {"etiqueta": "[3 2 1] — Dirección asimétrica (alta actividad de deslizamiento)", "uvw": [3, 2, 1]},
        {"etiqueta": "[2 2 3] — Dirección oblicua de alto índice", "uvw": [2, 2, 3]},
    ]

    # Precalcular ejemplos clave para validación instantánea
    ejemplos_schmid = {
        "FCC_321": calcular_schmid_direccion([3, 2, 1], sistemas_fcc)[:4],
        "BCC_321": calcular_schmid_direccion([3, 2, 1], sistemas_bcc)[:4],
        "FCC_100": calcular_schmid_direccion([1, 0, 0], sistemas_fcc)[:4],
    }

    # 3. Parámetros de difusión atómica (Ecuación de Arrhenius: D = D0 * exp(-Qd / (R * T)))
    mecanismos_difusion = [
        {
            "id": "C_en_Fe_gamma",
            "nombre": "Carbono intersticial en Fe-γ (Austenita FCC)",
            "descripcion": "Proceso industrial estándar de cementación (850–980 °C). Aunque la red FCC es más compacta (FEA=0.74), disuelve hasta 2.14% C y evita crecimiento excesivo de grano.",
            "D0": 2.3e-5,       # m^2/s
            "Qd": 148000.0,     # J/mol
            "T_min_c": 850,
            "T_max_c": 1000,
            "T_default_c": 925,
            "C0_default": 0.15, # %C núcleo (acero SAE 1015/8620)
            "Cs_default": 1.00, # %C atmósfera rica en carbono
        },
        {
            "id": "C_en_Fe_alpha",
            "nombre": "Carbono intersticial en Fe-α (Ferrita BCC)",
            "descripcion": "Difusión en red BCC más abierta (FEA=0.68): requiere menor energía de activación (80 kJ/mol), por lo que el carbono salta mucho más rápido, pero su solubilidad máxima es muy baja (0.022% C).",
            "D0": 6.2e-7,       # m^2/s
            "Qd": 80000.0,      # J/mol
            "T_min_c": 500,
            "T_max_c": 720,
            "T_default_c": 650,
            "C0_default": 0.005,
            "Cs_default": 0.022,
        },
        {
            "id": "N_en_Fe_alpha",
            "nombre": "Nitruración en Fe-α (Nitrógeno en Ferrita BCC)",
            "descripcion": "Endurecimiento superficial a baja temperatura (500–580 °C) sin necesidad de temple posterior, evitando distorsiones geométricas en engranes de precisión.",
            "D0": 3.0e-7,       # m^2/s
            "Qd": 76150.0,      # J/mol
            "T_min_c": 480,
            "T_max_c": 580,
            "T_default_c": 530,
            "C0_default": 0.01,
            "Cs_default": 0.60,
        },
    ]

    # 4. Tabla exacta de función error erf(z) para z in [0, 3.0] (paso 0.05) para validación en JS
    tabla_erf = [
        {"z": round(z, 2), "erf": round(math.erf(z), 6)}
        for z in np.arange(0.0, 3.05, 0.05)
    ]

    # 5. Retos interactivos para exponer en vivo frente al grupo (Pestaña 4)
    retos_clase = [
        {
            "id": "reto_cristal",
            "titulo": "Reto 1: Monocristal bajo Tensión Uniaxial",
            "pregunta": "Se aplica un esfuerzo axial de σ = 80 MPa en la dirección [3 2 1] a un monocristal con esfuerzo cortante crítico τ_CRSS = 25 MPa. Si el cristal es FCC ({111}<110>), ¿ocurre deslizamiento plástico? ¿Y si es BCC ({110}<111>)?",
            "pista": "Compara τ_R = σ · cos(φ) · cos(λ) del sistema más favorable contra τ_CRSS = 25 MPa.",
            "respuesta_corta": "En FCC SÍ desliza (τ_R = 27.99 MPa > 25 MPa en el sistema (111)[10-1], con m = 0.3499). En BCC NO desliza todavía (τ_R máximo = 23.33 MPa < 25 MPa en (110)[-11-1], con m = 0.2916).",
        },
        {
            "id": "reto_difusion",
            "titulo": "Reto 2: Ingeniería de Cementación en Engranes",
            "pregunta": "Un engrane de acero SAE 1015 (C₀ = 0.15%) entra a un horno de cementación con atmósfera Cs = 1.00% C durante t = 6 horas. Si queremos que a una profundidad de x = 0.60 mm alcance al menos Cx = 0.40% C para poder templarlo, ¿basta con hornear a 880 °C o necesitamos subir a 930 °C?",
            "pista": "Usa DiffusioLab (Pestaña 3) con C₀ = 0.15%, Cs = 1.00%, t = 6 h y prueba ambas temperaturas.",
            "respuesta_corta": "A 880 °C la capa efectiva (0.40% C) apenas llega a ~0.46 mm (insuficiente). Subiendo a 930 °C el coeficiente D se duplica por Arrhenius y la capa alcanza ~0.64 mm, cumpliendo el diseño.",
        },
    ]

    # 6. Metodología y fundamentos teóricos explicados de forma simple
    metodologia = {
        "titulo": "Fundamentos Físicos y Metodología de Cálculo",
        "resumen": "AtomForge modela dos fenómenos gobernados por defectos cristalinos (Unidades 3 y 4): el deslizamiento de dislocaciones mediante la Ley de Schmid y la difusión atómica asistida térmicamente mediante Arrhenius y la Segunda Ley de Fick.",
        "pasos": [
            {
                "paso": "1. Generación de Familias Cristalográficas",
                "descripcion": "A partir de los índices de Miller se construyeron computacionalmente los 12 sistemas compactos de deslizamiento para redes FCC ({111}<110>) y los 12 sistemas primarios para BCC ({110}<111>), verificando ortogonalidad n · s = 0.",
            },
            {
                "paso": "2. Ley de Schmid y Esfuerzo Cortante Resuelto",
                "descripcion": "Para cualquier eje de carga [u v w], se calculan los cosenos directores cos(φ) respecto a la normal del plano y cos(λ) respecto a la dirección de deslizamiento. El sistema activo es el que maximiza m = cos(φ)cos(λ) y supera τ_CRSS.",
            },
            {
                "paso": "3. Coeficiente de Difusión por Arrhenius",
                "descripcion": "La movilidad atómica por vacancias e intersticios se calcula con D(T) = D₀ · exp(-Qd / (R · T)), usando constantes experimentales de carbono y nitrógeno en hierro FCC y BCC (R = 8.314 J/mol·K).",
            },
            {
                "paso": "4. Perfil de Concentración (Segunda Ley de Fick)",
                "descripcion": "Bajo condición de fuente superficial constante Cs en un sólido semi-infinito, la concentración a profundidad x y tiempo t se resuelve exactamente mediante (Cs - C(x,t))/(Cs - C₀) = erf(x / (2√(D·t))).",
            },
        ],
    }

    payload = {
        "autores": ["Basurto Valdivieso Erick", "Coronel Reveles Guadalupe"],
        "constante_R": 8.314,
        "metales_base": metales_base,
        "sistemas_cristalinos": {
            "FCC": sistemas_fcc,
            "BCC": sistemas_bcc,
        },
        "direcciones_preset": direcciones_preset,
        "ejemplos_schmid": ejemplos_schmid,
        "mecanismos_difusion": mecanismos_difusion,
        "tabla_erf": tabla_erf,
        "retos_clase": retos_clase,
        "metodologia": metodologia,
    }

    with open(salida, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    print(f"Listo: '{salida.name}' generado en {carpeta}")
    print(f"- Sistemas FCC: {len(sistemas_fcc)} | Sistemas BCC: {len(sistemas_bcc)}")
    print(f"- Ejemplo [321] FCC max m: {ejemplos_schmid['FCC_321'][0]['m']}")


if __name__ == "__main__":
    main()