# Pipe Pressure Loss Uncertainty GUI (MEE Fluids Lab 2026)

Desktop GUI application for the **MEE Fluids Lab** to explore how measurement uncertainty assumptions affect the agreement between:

- **Experimental Darcy–Weisbach friction factor** (from measured pressure loss and flow), and
- **Theoretical Darcy–Weisbach friction factor** (from Colebrook–White).

Students paste tabular data copied from Excel, adjust uncertainty assumptions, and re-run the analysis to see when the **experimental uncertainty bounds** cover the **theoretical uncertainty bounds**.

---

## What this app is for

This tool supports a lab activity where students investigate the effect of measurement uncertainty assumptions on the estimated friction factor in pipe flow.

Students can change only the following uncertainty inputs:

- **Pressure loss uncertainty, ± (%)**  *(percentage)*
- **Volume uncertainty, ± (L)**  *(absolute uncertainty in litres)*
- **Time uncertainty, ± (s)**  *(absolute uncertainty in seconds)*

The app then performs Monte Carlo uncertainty propagation and visualizes:

- Experimental uncertainty band (e.g., p05–p95)
- Theoretical (Colebrook–White) uncertainty band
- Whether the theoretical interval is fully inside the experimental interval
- Whether the intervals overlap at all

---

## Features

- Desktop GUI (PySide6 / Qt for Python)
- Paste measurement rows directly from Excel (tab-separated) (Example measurements are provided in "Pressure loss in pipes - measurements.xlsx")
- Student-friendly sliders/spinboxes for uncertainty inputs
- Monte Carlo uncertainty propagation
- Experimental vs theoretical friction factor uncertainty plot
- Results summary and coverage metrics
- CSV export of results
- Windows executable can be distributed to students (no Python installation needed)

---

## Input format (paste from Excel)

Paste rows in this **exact tab-separated format**:

```text
D_mm    start_V_L    end_V_L    delta_t_s    delta_p_mmH2O
```

### Notes
- `D_mm` = pipe diameter in **mm**
- `start_V_L`, `end_V_L` = start/end volume readings in **L**
- `delta_t_s` = elapsed time in **s**
- `delta_p_mmH2O` = pressure loss in **mmH2O**

The app internally computes:

- `V_accumulated_L = end_V_L - start_V_L`

---

## Physics / calculation notes

### Reynolds number (implemented form)

Theoretical friction factor uses Reynolds number in the dynamic viscosity form:

\[
Re = \frac{\rho V D}{\mu}
\]

where:

- \( \rho \) = fluid density (kg/m³)
- \( V \) = velocity (m/s)
- \( D \) = pipe diameter (m)
- \( \mu \) = dynamic viscosity (Pa·s)

This is physically equivalent to:

\[
Re = \frac{V D}{\nu}
\]

because:

\[
\nu = \frac{\mu}{\rho}
\]

### Default fluid property values (editable in advanced settings)
- `rho = 1000 kg/m³`
- `mu = 0.0009764 Pa·s`

---

## Repository structure

```text
pipe-pressure-loss-gui/
├─ src/
│  └─ pipe_pressure_loss_gui.py
├─ README.md
├─ LICENSE
├─ environment.yml
├─ requirements.txt              # optional (pip users)
├─ requirements-dev.txt          # optional (pip + pyinstaller)
├─ CIV107 PressureLossinPipes.pdf (Lab document)
├─ Pressure loss in pipes - measurements.xlsx (example measurement data which can be used to copy and paste)
└─ .gitignore
```

---

## Option 1 — Use the prebuilt Windows executable (recommended for students)

If a prebuilt executable is provided in the **GitHub Releases** page:

1. Go to **Releases**
2. Download `PipePressureLossGUI.exe` (or the zipped `PipePressureLossGUI` folder if using `--onedir`)
3. Double-click to run

### No Python installation is required.

### If Windows blocks the file
Sometimes Windows marks downloaded executables as untrusted.

Try:
- Right-click the `.exe` → **Properties**
- Tick **Unblock** (if shown)
- Click **Apply**
- Run again

If your institution uses managed lab PCs, IT may need to whitelist the file/folder.

---

## Option 2 — Run from source (Conda environment)

This is for instructors/developers who want to run the Python source code directly.

### 1) Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/pipe-pressure-loss-gui.git
cd pipe-pressure-loss-gui
```

### 2) Create the Conda environment

```bash
conda env create -f environment.yml
```

### 3) Activate the environment

Use the environment name defined in `environment.yml` (example shown as `labgui`):

```bash
conda activate labgui
```

### 4) Run the app from source

```bash
python src\pipe_pressure_loss_gui.py
```

> On macOS/Linux, use `/` instead of `\`:
> `python src/pipe_pressure_loss_gui.py`

---

## Rebuild the Windows executable (PyInstaller)

This section is for users who want to repackage the app into a Windows `.exe`.

## Important
Build the Windows executable on **Windows**. PyInstaller is not a cross-compiler.

### 1) Activate the build environment

```bash
conda activate labgui
```

### 2) Install PyInstaller (if not already included)

If your `environment.yml` or `requirements-dev.txt` already includes PyInstaller, you can skip this.

```bash
pip install pyinstaller
```

### 3) Build a single-file executable (`--onefile`)

```bash
pyinstaller --noconfirm --onefile --windowed --name PipePressureLossGUI src\pipe_pressure_loss_gui.py
```

Output:

- `dist\PipePressureLossGUI.exe`

### 4) Alternative build (recommended for Qt apps): `--onedir`

Qt apps sometimes start faster / behave more reliably with `--onedir`.

```bash
pyinstaller --noconfirm --onedir --windowed --name PipePressureLossGUI src\pipe_pressure_loss_gui.py
```

Output:

- `dist\PipePressureLossGUI\`

### Which should I use?
- **Students / simple distribution:** `--onefile`
- **Best reliability / faster startup:** `--onedir`

---

## Example Windows rebuild workflow (from scratch)

```bash
conda env create -f environment.yml
conda activate labgui
python src\pipe_pressure_loss_gui.py
pyinstaller --noconfirm --onefile --windowed --name PipePressureLossGUI src\pipe_pressure_loss_gui.py
```

---

## Dependency management

This project can be reproduced with:

- **Conda**: `environment.yml` (recommended for rebuilders)
- **pip**: `requirements.txt` / `requirements-dev.txt` (optional convenience)

If using pip instead of conda:

```bash
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python src\pipe_pressure_loss_gui.py
```

For packaging builds:

```bash
pip install -r requirements-dev.txt
```

---

## Typical classroom deployment workflow (recommended)

1. Build the executable on a Windows machine
2. Test it on one **clean** lab PC
3. Copy the executable (or `--onedir` folder) to all lab PCs
4. (Optional) Create a desktop shortcut
5. Keep a backup copy on USB / network drive / GitHub Releases

This avoids needing Python, conda, or pip on student lab machines.

---

## Troubleshooting

### The app opens but plotting fails / looks broken
Make sure the executable was built from a tested environment and try `--onedir` packaging (often more robust for PySide6 + Matplotlib).

### The executable does not open on a lab PC
- Check Windows “Unblock” in file properties
- Move it to a local folder (e.g., `C:\FluidsLab\`)
- Ask IT to whitelist the app/folder if institution policy blocks unsigned executables

### Paste fails
Check that the pasted columns are in the exact order:

```text
D_mm    start_V_L    end_V_L    delta_t_s    delta_p_mmH2O
```

and that `delta_p` is in **mmH2O**, not mH2O.

---

## Author / contact

Developed by **Mustafa Onur Onen** using OpenAI ChatGPT v5.1  
moonen1@sheffield.ac.uk

For use in **MEE Fluids Lab 2026**.

---

## License

This project is licensed under the **MIT License**.

MIT License

Copyright (c) 2026 Mustafa Onur Onen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
