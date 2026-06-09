# Kaosa eta hurbilketa-unearen sentikortasuna

Dokumentu honek Apophis-en 2029ko Lur-hurbilketaren **sentikortasun kaotikoa** aztertzen du:
GrAL eta ASSIST-en arteko aldea hurbilpen aurretik **metro batzuetakoa** da, baina urtebetera
**milaka kilometrora** hazten da. Hemen frogatzen da hori ez dela ez zenbakizko akatsa ez eredu-akatsa,
baizik eta Lur-hurbilketaren **sentikortasun fisikoa**, eta hori **neurgarria** dela.

Erreferentziak: `scripts/error_budget.py`, `scripts/extract_assist_states.py`,
`scripts/multistart_encounter.jl`, `notebooks/GrAL_vs_ASSIST_Konparaketa.ipynb` (15. atala).

---

## 1. Tolerantzia eta zenbakizko zehaztasuna

- **Tolerantzia berdina** integratzaile bietan: IRKGL16 eta Vern9, `reltol = abstol = 1e-12`.
  Bi integratzaile alderatzeko, tolerantzia berdina ezinbestekoa da.
- Tolerantzia **urratseko errore lokala** kontrolatzen du, ez globala; errore globala pilatu egiten da
  eta hurbilketak anplifikatzen du.
- **Nahikoa da:** hurbilpen aurretik IRKGL16 ↔ Vern9 = **0,012 m** (τ=1e-12), eta τ=1e-13-ra estutzeak
  ozta-ozta hobetzen du (0,009 m — Float64-ren zorua). 0,012 m ≪ eredu-diferentziak (metroak–ehunka metro),
  beraz tolerantzia ez da muga.

---

## 2. Kaosa neur daiteke? Bai.

Hurbilketa-unean **posizio-perturbazio kontrolatu** bat sartu (along-track norabidean) eta
urte amaierara arteko hazkundea neurtu (GrAL+Yarkovsky, τ=1e-12):

| δ hurbilketan | Δ urte amaieran | irabazia (Δ/δ) |
|---|---|---|
| 1 m    | 4,2 km   | 4187× |
| 10 m   | 41,9 km  | 4189× |
| 100 m  | 418,9 km | 4189× |
| 1000 m | 4189 km  | 4189× |

Irabazia **konstantea** da (~4189×) → erantzun **lineala** → neurketa ondo definitua.

> **Apophis-en hurbilketak along-track posizio-diferentzia bat ~4200× anplifikatzen du**
> hurbilketatik urte amaierara. Hori da "kaosaren" neurri kuantitatiboa kasu honetan.

Esponentzial gisa hartuta, e-toleste denbora ~32 egunekoa litzateke, baina **ez da Lyapunov denbora
egokia** — ikus hurrengo atala.

---

## 3. Mekanismoa: hurbilketa anplifikadore, gero drift lineala

"Kaosa" ez da urte osoan zeharreko hazkunde esponentzial leuna; **Lur-hurbilketa bakar batek**
anplifikadore diskretu gisa jokatzen du:

- **Hurbilpen aurretik:** hazkunde leuna.
- **Hurbilketa-unea (2029-04-13):** ia berehalako **~4200× anplifikazioa**. Hemen dago sentikortasun
  kaotikoaren muina (b-planoaren geometria: huts-distantziaren aldaketa txikia → orbita-aldaketa handia).
- **Hurbilpen ostean:** hazkundea **LINEALA** da denboran (along-track drift, periodo-aldaketaren ondorioz),
  ez esponentziala. (Egiaztapena: 82→562→2102 km, 109→182→365 egunetan; faktorea 3,7 ≈ denbora-erlazioa 3,5.)

Beraz "kaosa" hurbilketa-unearen sentikortasunean dago kontzentratua; urteko anplifikazioa = hurbilketaren
irabazia × ondoko drift lineala.

---

## 4. Hurbilketa-uneko egoeratik abiatzea — eta harmoniken zeinu-akatsa aurkitzea

ASSIST eredu osoaren **egoera zehatza** (posizioa **eta** abiadura) hartu eta GrAL+Yarkovsky-rekin 2030 arte
integratu. Hasieran esperimentu honek urte-amaierako aldea **oso handia** zela erakutsi zuen (~3700 km), eta
horri jarraituz **harmoniken zeinu-akats** bat agerrarazi zen GrAL-en: J2/J3/J4 blokeek `d = gorputza − ast`
erabiltzen zuten, ASSIST-en `d = ast − gorputza` ordez; J2/J4 azelerazioa **d-rekiko bakoitia** denez, zeinua
alderantzizkatua zegoen (zenbakizko proba: ratio −1). Akatsa **konponduta** (`src/equations.jl`-eko lau funtziotan):

**Perigeoa:** t = 10695,408 egun, geozentriko distantzia = **38.015 km** (~6 Lur-erradio), 2029-04-13.

| Hasiera-unea | Δ 2030 (akatsarekin) | Δ 2030 (**konponduta**) | Perigeoa zeharkatzen? |
|---|---|---|---|
| Holman IC (100 egun lehenago) | 3699 km | **30,8 km** | Bai |
| Perigeoa − 1 egun | 3699 km | 33,5 km | Bai |
| Perigeoa | 1584 km | 8,3 km | (unean) |
| **Perigeoa + 1 egun** | 319 m | **106 m** | Ez |

**Emaitza nagusiak (konponduta):**

1. GrAL+Yarkovsky vs ASSIST ona = **30,8 km** urtebetera — ASSIST↔JPL zoruaren (**38,66 km**) parekoa.
2. **Yarkovsky-k aldea 1619 km → 30,8 km murrizten du**, ASSIST-en jokabide berdina (1611 → 38,66 km).
3. Perigeoa zeharkatzeak oraindik anplifikatzen du (106 m ostean → 30 km aurretik), baina orain hazia **txikia**
   da (integratzailea + GR_EIH eredu-aukera), **ez harmoniken akatsa**.

---

## 5. Ondorioak (tesirako)

- Hurbilketak edozein hazi **~4200× anplifikatzen du** (2.–3. atalak): hori **fisika erreala** da, aldagaiekiko
  independentea.
- **Baina** lehen ikusitako urte-amaierako milaka km-ko aldea **zati handian zenbakizko artefaktua** zen
  (harmoniken zeinu-akatsa), **ez kaos fisiko hutsa**. Konponduta, aldea **30 km**-ra jaisten da (ASSIST↔JPL parekoa).
  Hortaz, kaos-anplifikazioak hazi **txikia** behar du; akats sistematiko batek hazia handitzen badu, emaitza
  okerra anplifikatzen da.
- Integratzailea balioztatuta dago: IRKGL16 ↔ Vern9 = **3,7 m** urtebetera (indar berdinekin).
- **Yarkovsky orain urte-amaieran ere ageri da** (1619 → 30,8 km), ASSIST bezala. Hurbilpen aurretik
  GrAL+Yarkovsky-k eredu osoa **~2,7 m**-ra jarraitzen du 100. egunean.
- Geratzen den **~30 km** zorua GR_SIMPLE vs GR_EIH ereduaren + integratzaile-haziaren ondorioa da,
  hurbilketak anplifikatuta — ASSIST↔JPL-en mailakoa.
- GrAL-i Yarkovsky gehitu zaio (`physics.yarkovsky`), eta horrek eredu-hutsune nagusia ixten du
  hurbilpen aurretik (~400 m → ~2 m).
