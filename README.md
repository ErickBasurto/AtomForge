# AtomForge — Laboratorio Interactivo de Cristalografía y Difusión Atómica

Plataforma web interactiva orientada a visualizar y analizar cómo los defectos cristalinos gobiernan el comportamiento mecánico y térmico de los metales a escala atómica, integrando modos de visualización dual (**Modo Explorador** intuitivo y **Modo Ingeniero** técnico).

🔗 **Ver aplicación interactiva en vivo:** [https://erickbasurto.github.io/AtomForge/](https://erickbasurto.github.io/AtomForge/)

---

## Desarrolladores
* **Basurto Valdivieso Erick**
* **Coronel Reveles Guadalupe**

---

## Características Principales y Módulos

* **Selector Global de Audiencia (Modo Explorador / Modo Ingeniero):** Alterna en tiempo real entre explicaciones cotidianas con indicadores visuales tipo semáforo y el formalismo matemático completo (cosenos directores, constantes de Arrhenius y solución analítica de Fick).
* **1. Mapa Atómico:** Contraste entre defectos lineales (dislocaciones) y defectos puntuales (vacancias e intersticios), catálogo de metales base ($\text{Cu}$, $\text{Al}$, $\text{Ni}$, $\text{Fe-}\alpha$, $\text{Fe-}\gamma$) y casos de estudio aplicados en **ingeniería aeroespacial** (álabes monocristalinos en `[001]`) e **ingeniería biomédica** (prótesis ortopédicas de titanio y stents cardiovasculares).
* **2. CrystalSlip 3D (Ley de Schmid):** Visualizador 3D rotable de celdas unitarias **FCC** (`{111}<110>`) y **BCC** (`{110}<111>`) que calcula el factor de Schmid ($m = \cos\phi \cdot \cos\lambda$) para cualquier dirección de carga `[u v w]` y anima el deslizamiento atómico al superar $\tau_{\text{CRSS}}$.
* **3. DiffusioLab (Arrhenius y Segunda Ley de Fick):** Simulador de tratamientos termoquímicos (cementación y nitruración) que combina una animación atómica en `HTML5 Canvas` con la curva de concentración $C(x, t)$ y cálculo de profundidad efectiva de capa endurecida.
* **4. Laboratorio de Desafíos (Minijuegos):** Dos simulaciones interactivas de evaluación (*Misión: Orienta el Cristal* y *Maestro del Horno de Cementación*) con retroalimentación física inmediata.

---

## Flujo de Trabajo y Reproducibilidad

1. **Cálculo Cristalográfico y Termodinámico (`generar_datos_atomicos.py`):** Genera en Python los 24 sistemas de deslizamiento de FCC y BCC, las constantes de difusión ($D_0, Q_d$) y los valores de referencia de la función error $\text{erf}(z)$ en `datos_atomicos.json`.
2. **Ensamblaje Web (`armar_web.py` $\rightarrow$ `index.html`):** Inyecta `datos_atomicos.json` en la interfaz construida con `HTML5`, `Tailwind CSS`, `Plotly.js` y `Canvas`, produciendo un único archivo `index.html` autocontenido.

---

## Estructura del Repositorio
* `index.html` — Aplicación web interactiva desplegada en GitHub Pages.
* `generar_datos_atomicos.py` — Script en Python (`numpy`, `math`) para cálculos cristalográficos y de difusión.
* `datos_atomicos.json` — Constantes físicas, sistemas cristalinos y parámetros exportados.
* `armar_web.py` — Script de ensamblaje automático entre datos y plantilla web.
