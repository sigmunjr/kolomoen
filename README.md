# Laveste terrengpassasje: Heradsbygd–Kolomoen

Pakken inneholder ferdig beregnede resultater, et selvstendig Python-skript og en
installerbar skrivebordsviser bygget med Node.js/Electron, og en statisk nettviser
som kan publiseres på GitHub Pages.

**Resultat for den opplastede filen:** minste sammenhengende rasterterskel **222,5 m**;
representativ korteste rute ved denne terskelen **29,106 km**. Ved **222,4 m** er
endepunktene adskilt. Det avgjørende rasterområdet ligger rundt **60.79122 N, 11.52631 E**.
Se `results/RESULTAT.md` for metode, metadata, høydeforbehold og hva dette ikke beviser.

## GitHub Pages

Nettviseren publiseres på **https://sigmunjr.github.io/kolomoen/**.
Hele 3D-viseren, egne terrengdata, nedlastinger og alle tre høydeprofiler fungerer
direkte i nettleseren. Det krever verken installasjon, Python-server eller
Cesium ion-token. ArcGIS og OpenStreetMap hentes fortsatt fra sine nettjenester.

### Første gangs oppsett

1. Push prosjektet og `.github/workflows/pages.yml` til `main` på GitHub.
2. Åpne [Settings → Pages](https://github.com/sigmunjr/kolomoen/settings/pages).
   Under **Build and deployment → Source**, velg **GitHub Actions**.
3. Åpne [Actions](https://github.com/sigmunjr/kolomoen/actions), velg
   **Deploy viewer to GitHub Pages**, og kjør **Run workflow** dersom en kjøring
   ikke allerede er i gang. Vent til både `build` og `deploy` er grønne.
4. Åpne https://sigmunjr.github.io/kolomoen/.

Ingen egne secrets, tilgangstokener, ion-konto eller `gh-pages`-gren trengs.
Workflowen bruker GitHubs innebygde `GITHUB_TOKEN`. Senere push til `main`
bygger og publiserer nettstedet automatisk. Pull requests bygger og validerer
nettstedet uten å publisere det.

### Bygg og prøv lokalt

```bash
npm ci --ignore-scripts
npm run build:pages
python3 -m http.server 8000 --bind 127.0.0.1 --directory dist
```

Åpne http://localhost:8000/pages/. `dist/pages` inneholder den ferdige statiske
siden. Alle adresser, også Cesiums worker- og teksturfiler, er relative til siden
og fungerer under prosjektstien `/kolomoen/` eller på et eget domene.
Bygget inkluderer Cesiums lisens og de nødvendige `Workers`, `Assets`,
`ThirdParty` og `Widgets`, men utelater den eldre innebygde HTML-filen og
analyse-GeoTIFF-en. Terrengdataene må være faktiske filer, ikke Git LFS-pekere;
workflowen henter LFS-data dersom det senere tas i bruk.

Ved feil: se `build`/`deploy`-loggen i Actions og kontroller at Source fremdeles
er **GitHub Actions**. Hvis kartet er tomt, bruk en nettleser med WebGL og
maskinvareakselerasjon. Bakgrunnskart kan slås av for å bruke bare egne terrengdata.

## 1. Installer eller start skrivebordsappen

På Linux kan den ferdig bygde `dist/Heradsbygd Kolomoen-1.1.0-linux-x86_64.AppImage`
startes direkte uten Node.js; gjør filen kjørbar dersom filbehandleren ber om det.
Appen pakker CesiumJS og alle terreng- og analysedata lokalt. Bare de valgfrie
bakgrunnskartene fra ArcGIS eller OpenStreetMap krever internettilgang.

For å starte fra kildekoden må Node.js og npm være installert:

```bash
npm install
npm start
```

Lag en installasjonspakke for den aktuelle plattformen med én av disse kommandoene:

```bash
npm run dist:linux   # AppImage
npm run dist:win     # NSIS-installasjonsprogram
npm run dist:mac     # DMG
```

Windows- og macOS-pakkene bør bygges på henholdsvis Windows og macOS. Bygget kopierer
den versjonslåste CesiumJS-distribusjonen til appen automatisk.

Python-serveren er fortsatt et lett alternativ for utvikling:

Pakk ut ZIP-filen, åpne en terminal i `heradsbygd_kolomoen`, og kjør:

```bash
python3 open_viewer.py
```

Nettleseren åpnes på `http://localhost:8000`. På Windows kan kommandoen hete
`python open_viewer.py` eller `py open_viewer.py`.

Alternativt, uten å åpne nettleseren automatisk:

```bash
python3 -m http.server 8000 --bind 127.0.0.1 --directory results
```

Åpne så `http://localhost:8000` manuelt. Den oppdelte `index.html` må få lese
terrengfilen og JSON-filene over HTTP. Serveren eksponeres bare på din egen maskin.
Avslutt med Ctrl+C. En opptatt port kan byttes med `python3 open_viewer.py --port 8001`.

Selve analysen, CesiumJS **1.145** og alle terrengdata er lokale. **Ingen Cesium
ion-konto eller token brukes.** ArcGIS World Imagery og OpenStreetMap kan velges
som bakgrunnskart og bruker nett. Ingen høydedata eller rutefiler lastes opp til
Cesium ion.

Den tidligere selvstendige HTML-filen beholdes bare som eldre eksport. Bruk
skrivebordsappen eller den lokale Python-serveren for pålitelig lasting av de store
terrengfilene og Cesiums `Workers`, `ThirdParty`, `Assets` og `Widgets`.

### Kontroller i kartet

Velg terskelnivået under, ved eller over 222,5 m. Fargeflatene viser sammenhengende
terrengceller under nivået, ikke et estimat av faktisk oversvømmelse. Slå av/på
korteste rute, rett referanse og dalvektet alternativ. Bruk «Til terskelen» for
nærbilde. Høydeforstørrelse 1×, 3× eller 6× påvirker bare visningen.
Profilen nederst bruker originale ruteceller. Velg mellom laveste passasje, rett
referanselinje og **Dal nær (Indigo Valley Line)**; pek på profilen for å markere stedet.
Bakgrunnskartet kan slås av og på, byttes mellom ArcGIS World Imagery og
OpenStreetMap og blandes med den lokale terrengskyggen med glidebryteren.

## 2. Gjenta analysen

Originalfilen `6706_3_10m_z32.tif` er **ikke duplisert i ZIP-en**. Legg TIFF-en ved
siden av `lowest_passage.py`, eller bruk en full filsti. Originalen endres ikke.

Linux/macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python lowest_passage.py 6706_3_10m_z32.tif --out results
python open_viewer.py
```

Windows PowerShell:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe lowest_passage.py 6706_3_10m_z32.tif --out results
.\.venv\Scripts\python.exe open_viewer.py
```

Utdata i `--out` overskrives ved ny kjøring. Bruk et nytt mappenavn for å bevare en
eldre beregning. `lowest_passage.py` inneholder også hele HTML/JavaScript-malen;
filen kan brukes selvstendig etter installasjon av bibliotekene.

Koordinatene fra brukerens melding er standardverdiene:

```bash
python lowest_passage.py 6706_3_10m_z32.tif \
  --start 60.820232979269655 11.635771799665232 \
  --end 60.69954161460697 11.265232590579124 \
  --out results
```

Argumentrekkefølgen er alltid **LAT LON**. GeoJSON bruker standardrekkefølgen
**LON LAT**. Skriptet håndterer denne forskjellen eksplisitt.

### Andre avgrensninger og utvalg

`--via LAT LON` kan gjentas for å kreve at ruten passerer bestemte, eksplisitte
celler, for eksempel et kvalitetssikret punkt i et bestemt passasjebelte ved
Rokosjøen. Ingen slik tvang er lagt på standardberegningen. Med flere punkter
beregnes først en felles global terskel; deretter korteste rute via alle punktene
under den samme terskelen.

`--search-buffer-km 5` begrenser analysen til et rektangel som omfatter alle
oppgitte punkter med fem kilometers margin. Dette endrer søkedomenet og kan
endre terskelen; en avgrenset kjøring kan ikke brukes som bevis for global optimalitet.
Standard er **hele den opplastede rasterflisen**.

`--valley-weight 0` deaktiverer det illustrative dalvektede alternativet.
`--profile-step-m 100` velger 100 m mellom profilprøvene og gir filnavnene
`route_100m.csv` og `reference_100m.csv` i stedet for `50m`-filene.
`--skip-web` utelater web-eksport. `--web-max-side 3072` gir et tettere visningsraster,
ikke en bedre vitenskapelig analyse.

25,4 millioner celler i denne filen krever omtrent 1–2 GB ledig RAM, avhengig av
miljø og implementasjonsdetaljer. Første kjøring kompilerer A*-kjernen med Numba.
1-metersdata over like stort areal er **100 ganger så mange celler** og må behandles
med en passende mindre analyseavgrensning, mer RAM eller en separat flislagt løsning.
Skriptet nekter som standard å lese mer enn 50 millioner analyseceller og vil ikke
uten beskjed redusere oppløsningen. `--max-cells` kan endres ved tilstrekkelig minne.

## 3. Filer

| Fil | Innhold |
|---|---|
| `lowest_passage.py` | Hele analysen, profil-/karteksport og innebygd Cesium-malside |
| `desktop/main.js`, `package.json` | Electron-app, lokal filserver og installasjonsbygg |
| `open_viewer.py` | Alternativ lokal HTTP-server og oppstart av nettleser; kun standardbibliotek |
| `results/RESULTAT.md` | Resultat, metode og begrensninger |
| `results/summary.json` | Filhash, CRS, rasteroppløsning, koordinater og validering |
| `results/index.html` | CesiumJS-side som leser de lokale resultatfilene |
| `results/Heradsbygd-Kolomoen.html` | Eldre selvstendig HTML-eksport; skrivebordsappen anbefales |
| `results/route_all_cells.csv` | Alle 2 569 ruteceller med UTM, lat/lon, høyde og avstand |
| `results/route_50m.csv` | Prøvepunkter hver 50. meter langs den valgte ruten |
| `results/reference_27_NEW_points.csv` | 27 **nye** punkter langs den rette referanselinjen, med rasterhøyder |
| `results/reference_50m.csv` | 487 punkter langs rett referanselinje |
| `results/reference_all_crossed_cells.csv` | Kontroll mot hver celle krysset av den rette linjen |
| `results/critical_plateaus.csv` | Kritisk terskelområde og andre høydekoteberøringer, eksplisitt skilt |
| `results/critical_plateaus.geojson` | Geometrien til disse rasterområdene |
| `results/connected_at_threshold.tif` | Sammenhengende lavområde ved terskelen, i originaloppløsning |
| `results/connectivity_checks.csv` | Sammenheng før, ved og over terskelen |
| `results/saddle_north_south_2km.csv` | Fast nord–sør-profil gjennom representativ terskelcelle |
| `results/routes.geojson` | Ruter og punkter, 2D-geometri; høyder i egenskapene |
| `results/overview.png` | Terrengkart med ruter og terskelområde |
| `results/profiles.png` | Profil langs optimal rute sammenlignet med rett linje |
| `results/saddle_section.png` | Tverrprofil gjennom terskelområdet |
| `results/terrain.f32`, `hillshade.png`, `connected_*.png`, `viewer_data.json` | Visningsdata; ikke erstatning for originalrasteret |

GeoJSON-geometrien er med hensikt 2D: tredje koordinat i vanlig GeoJSON ville
blitt tolket som en ellipsoidehøyde. De utransformerte høydene er derfor lagret
under `native_elevations_m` i stedet. CSV bruker meter i originalrasterets
høydereferanse. CSV-rad/kolonne er nullbasert i originalfilen.

Den innebygde WKT-en viser en uoverensstemmelse mellom EPSG-ID 3044 og den oppgitte
øst/nord-akseordenen. Den faktiske projeksjonen/akseordenen identifiseres av pyproj
som EPSG:25832 (UTM 32). Både original-WKT og avledet identifikasjon bevares i metadata.
Skriptet bruker `always_xy=True` og behandler ikke dataene som UTM 33.

## 4. Vitenskapelige forbehold

PDF-en beskriver behov for en sammenhengende passasje, terskler og profiler som
forarbeid til en geologisk vurdering. Denne pakken beregner topografiske produkter,
men gjør ikke hele den historiske/hydrologiske rekonstruksjonen.

Rasterfilen er **10 m**, ikke 1 m. PDF-ens omtrentlige endepunkter og historiske
27 punkter er ikke brukt som om de var identiske med de nye, presise koordinatene.
De 27 CSV-punktene er et nytt, eksplisitt definert datasett.

Analysen bruker rå høydedata uten senkefylling, elveinnbrenning, terrengglatting eller
manuell fjerning av vegfyllinger. Ved minimumsterskelen kan passasjen ha null
effektiv vanndybde. Et faktisk avløp avhenger blant annet av vannstand, tilgjengelig
vann, fall i vannspeilet, andre utløp og eventuelle is-/terrengbarrierer.

Begge endepunktenes lavområder når yttergrensen til rasteret allerede under
terskelen. Derfor er en lavere omvei **utenfor** filens dekningsområde ikke utelukket.
Ruten har lavest maksimal rasterhøyde innenfor den analyserte flisen, ikke bevist
lavest terskel på et større, ubegrenset kart.

Det vertikale datumet er ikke oppgitt entydig i TIFF-en. Absolutt geodetisk
plassering i Cesium krever derfor en senere bekreftet høydetransformasjon.
Viseren bruker en konsistent **lokal visningsreferanse**, ikke en påstått korrekt
NN2000→ellipsoidekonvertering. Rutehøydene endres ikke av hvordan kartet vises.

## 5. Tester og programstatus

```bash
python lowest_passage.py --self-test
```

Algoritmen er kjørt på den opplastede filen. Testene omfatter 120 tilfeldige små
rastere sammenlignet med uttømmende terskelkontroll, NoData-barrierer, negative
høyder, diagonal hjørnekontakt, endepunktdominerte terskler, korteste rute og via-punkter.
Både 4-naboskap og ukritisk 8-naboskap gir 222,5 m i denne filen.

Cesium-koden er syntakskontrollert og bygger på den offisielle API-dokumentasjonen.
En full nettlesertest mot den lokale webserveren ble **ikke fullført**: testmiljøets
nettleser blokkerte localhost med `ERR_BLOCKED_BY_ADMINISTRATOR`. Nettleserdelen må
derfor fortsatt funksjonstestes i et normalt lokalt nettlesermiljø. Dette er en
begrensning i testing av viseren, ikke i beregningen av terskelen eller CSV-/rasterdataene.

## Tekniske referanser

- Cesium CustomHeightmapTerrainProvider: https://cesium.com/learn/cesiumjs/ref-doc/CustomHeightmapTerrainProvider.html
- Cesium SingleTileImageryProvider: https://cesium.com/learn/cesiumjs/ref-doc/SingleTileImageryProvider.html
- Cesium Cartesian3 (høyde over ellipsoiden): https://cesium.com/learn/cesiumjs/ref-doc/Cartesian3.html
- CesiumJS-oppsett og lokale assets: https://cesium.com/learn/cesiumjs-learn/cesiumjs-quickstart/
- Rasterio reprojeksjon: https://rasterio.readthedocs.io/en/stable/topics/reproject.html

Disse referansene beskriver teknologien. De er ikke kilder til de beregnede
lokale høydene, som kommer fra den opplastede TIFF-en.
