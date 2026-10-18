# Heradsbygd–Kolomoen: laveste sammenhengende terrengpassasje

## Beregnet resultat

**Minste forbindende rasterhøyde: 222.5 m.**
Korteste rute under dette taket: **29.106 km**.
Ved forrige lavere rasterhøyde (222.4 m) er endepunktene ikke forbundet.
Dette er beregnet fra alle relevante celler i originalrasteret, ikke fra et anslag eller en 50-metersprofil.

| Størrelse | Verdi |
|---|---:|
| Startcellens høyde | 172.9 m |
| Sluttcellens høyde | 184.2 m |
| Laveste forbindende terskel | 222.5 m |
| Korteste rute ved denne terskelen | 29.106 km |
| Rett referanselinje, lengde | 24.263 km |
| Høyeste kryssede rastercelle langs rett linje | 345.1 m |

Den kritiske sammenkoblingen skjer i et rasterplatå med 46 celler på 222.5 m. Representativ celle: 60.79122374 N, 11.52630782 E; øst 637510.0, nord 6742180.0. Rasterstøtten er 110 × 120 m i utstrekning. Dette er ikke en målt kanalbredde eller en endelig geologisk sadel.

## Kilde og referansesystem

Fil: `6706_3_10m_z32.tif`. SHA256: `8175b418cbfbb758b635c0a85424637c0ee9d51af6c1c3d668859cb54e4d1e2a`.
Originalraster: 5041 × 5041 celler.
Cellestørrelse: 10 × 10 m.
CRS-navn: ETRS89 / UTM zone 32N (N-E). Innebygd CRS-ID: {'authority': 'EPSG', 'code': 3044}.
Pyproj identifiserer den faktiske projeksjonen/akseordenen som ('EPSG', '25832').
Rasterio returnerer ingen entydig autoritetskode for den opprinnelige WKT-en. Dette avviket
bevares i metadata i stedet for å skjules. Rasterets X/Y behandles som øst/nord;
koordinattransformasjoner bruker `always_xy=True`. `AREA_OR_POINT` er `Point`.
Rasterio/GDAL håndterer dette allerede: ingen ekstra halv-celle-forskyvning legges til.

**Vertikalt datum er ikke bekreftet av filens metadata.** Tallene er originalrasterets høyder i meter.
De må ikke automatisk omtales som kvalitetssikrede NN2000- eller ellipsoidehøyder.
Desimalen i en rasterverdi er ikke en påstand om tilsvarende målenøyaktighet.

## Definisjon og algoritme

For en sammenhengende rute P minimeres `max(z for cell in P)`.
Vi øker et høydenivå og finner det første nivået som kobler sammen endepunktene
(og eventuelle eksplisitte `--via`-punkter). Deretter søker A* den korteste ruten
innenfor dette nivået. Denne totrinnsmetoden unngår å blande høyde og lengde i ett
vilkårlig hovedmål. Alle sluttverdier valideres med sammenheng like under og ved terskelen.

Sammenheng beregnes med felles celleside. Ruten kan også gå diagonalt når **begge**
sidestøttende celler ligger under høydegrensen. Det hindrer forbindelser som bare
berører hverandre i et hjørne. Som følsomhetskontroll ga også ukritisk 8-naboskap
222.5 m. Dette er et resultat på rastergrafen,
ikke en eksakt kontinuerlig optimal rute over et ukjent terreng mellom målepunktene.

Det dalvektede alternativet bruker samme terskel og en sekundær kostnad
`ds * (1 + weight * ((max(z1,z2) - raster_min)/(ceiling-raster_min))**4)`.
Dette er en illustrativ preferanse for lavere terreng, ikke en hydraulisk strømmodell.
Ruter er ikke glattet eller forenklet før analyse/eksport.

## Avgrensning

Both subcritical endpoint components touch the analysis boundary or NoData. A lower route outside the supplied/selected domain is not ruled out.

Søkedomenet er 5041 × 5041 celler.
Ingen lavere passasje **innenfor dette domenet** er funnet eller tillatt av rastertesten.
Det er ikke satt et krav om passasje gjennom Rokosjøen, siden ingen koordinat eller
avgrenset passasjesone der ble oppgitt. `--via LAT LON` kan legge inn et slikt eksplisitt krav.

## Hva dette sier – og ikke sier – om vann

Terskelen er høyere enn begge endepunktene. Dermed finnes ingen hele veien
ikke-stigende terrengbunn mellom akkurat disse punktene i denne rastermodellen.
Det utelukker ikke vanntransport over en terskel ved høyere vannstand eller oppdemming.
Minimax-testen er symmetrisk og avgjør ikke hvilken retning vannet ville beveget seg.

Nivålagene i kartet er sammenhengende **terrengceller under et nivå**, ikke dokumenterte
innsjøer eller beregnede flomsoner. Ved selve minimumsterskelen kan forbindelsen ha
null vanndybde; vannføring krever et reelt tverrsnitt, hydraulisk fall og tilstrekkelig vann.
Tidligere isdammer, innløp, andre utløp, senere erosjon og avsetninger er ikke modellert.
Vannoverflater i DTM er ikke innsjøbunn. Vegfyllinger og små naturlige terskler kan ha
avgjørende betydning og må kontrolleres i bedre data og eventuelt i felt.

## Referansepunkter og profiler

`reference_27_NEW_points.csv` inneholder **27 nye punkter** jevnt fordelt langs en
rett UTM-linje mellom koordinatene i brukerens melding. Disse er ikke de historiske
27 punktene som PDF-en sier må kvalitetssikres. PDF-ens omtrentlige endepunkter
avviker også fra brukerens nye endepunkter. Se `input_points.csv` for faktiske innparametre.

`reference_50m.csv` og `route_50m.csv` gir profiler med valgt intervall (standard 50 m).
Høydene er hentet fra cellen ved prøvepunktet; de er ikke nye feltmålinger.
`reference_all_crossed_cells.csv` og `route_all_cells.csv` bevarer tettere rasterkontroll.
Høyeste verdi på referanselinjen er kontrollert mot kryssede celler, ikke bare 50 m-prøver.
Rutens lengde måles mellom cellesentre. Endepunktene kobles til sin inneholdende celle;
små forskyvninger er eksplisitt oppgitt som `snap_distance_m` i `input_points.csv`.

`critical_plateaus.csv` skiller det forbindende terskelområdet fra andre steder der
ruten bare berører samme høydekote. `saddle_north_south_2km.csv` er en fast nord–sør-
profil gjennom en representativ terskelcelle, ikke en ferdig tolkning av landskapsformen.

## CesiumJS

Kjør `python -m http.server 8000 --bind 127.0.0.1 --directory results` fra prosjektmappen,
og åpne `http://localhost:8000`.
Se `README.md` i prosjektmappen for oppsett.

Viseren bruker egne terrengdata fra TIFF-en, ikke Cesium World Terrain.
Analyse og profiler bruker originaloppløsningen; nettleserterreng er reprojisert og
forenklet kun for visning. Terreng, ruter og terskelmarkører har konsistent lokal
visningshøyde, men originalhøydene er **ikke** transformert til ellipsoidehøyder.
En korrekt geodetisk overlagring krever bekreftet vertikalt datum og riktig HREF/geoid-
transformasjon; en vilkårlig konstant forskyvning er ikke lagt inn som om den var korrekt.

Kartet kan vise en valgt, felles høydeforstørrelse; avlest analysehøyde er uforandret.
Utenfor det eksporterte utsnittet er det ingen terrengdata (visningsellipsoide ved 0).
GeoJSON-geometri er 2D for ikke å feilmerke rå rasterhøyder som WGS84-ellipsoidehøyder.

## Reproduserbarhet

Programversjon 1.0.0; Python 3.13.5.
`summary.json` inneholder metadata, filhash, valideringsresultater, koordinater og forbehold.
Kjør `python lowest_passage.py --self-test` for syntetiske enhetstester.
