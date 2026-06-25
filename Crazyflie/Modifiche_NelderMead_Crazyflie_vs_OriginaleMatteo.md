# Modifiche apportate per adattare la versione Crazyflie (STM32) della libreria Nelder-Mead rispetto all'originale di Matteo

Questo documento elenca le principali differenze e modifiche tra la versione originale della libreria Nelder-Mead di Matteo e quella adattata per l'utilizzo su Crazyflie (STM32).

## 1. Gestione della memoria
- **Crazyflie:** Uso di array statici e allocazione fissa per evitare allocazioni dinamiche (malloc/free) non raccomandate su microcontrollori.
- **Originale:** Possibile uso di allocazione dinamica o array di dimensione variabile.

## 2. Dipendenze e inclusioni
- **Crazyflie:** Inclusione di header specifici per STM32 e firmware Crazyflie, rimozione di dipendenze non compatibili.
- **Originale:** Inclusione di librerie standard C e header generici.

## 3. Tipi di dato
- **Crazyflie:** Uso esplicito di tipi come `float` o `double` in base alle limitazioni hardware.
- **Originale:** Possibile uso di tipi generici o configurabili.

## 4. Funzioni di utilità
- **Crazyflie:** Funzioni matematiche semplificate o riscritte per evitare funzioni complesse non ottimizzate per embedded (es. `pow`, `sqrt`).
- **Originale:** Uso diretto di funzioni matematiche standard.

## 5. Debug e logging
- **Crazyflie:** Rimozione di stampe su console (`printf`) o sostituzione con macro di debug compatibili con il firmware.
- **Originale:** Uso di `printf` per debug e log.

## 6. Parametri e costanti
- **Crazyflie:** Parametri e costanti definiti come macro o variabili globali per ottimizzare la gestione della memoria.
- **Originale:** Parametri definiti localmente o come variabili configurabili.

## 7. Strutture dati
- **Crazyflie:** Strutture dati semplificate e ottimizzate per l'uso su microcontrollore.
- **Originale:** Strutture dati più flessibili e generiche.

## 8. Gestione degli errori
- **Crazyflie:** Gestione degli errori semplificata, spesso tramite return code o macro.
- **Originale:** Gestione degli errori più dettagliata.

## 9. Ottimizzazioni specifiche
- **Crazyflie:** Ottimizzazioni per ridurre il consumo di memoria e CPU.
- **Originale:** Algoritmo più generico e meno ottimizzato per risorse limitate.

## 10. Interfaccia con il resto del firmware
- **Crazyflie:** Funzioni e variabili adattate per integrarsi con il firmware Crazyflie.
- **Originale:** Interfaccia generica, non pensata per firmware specifici.

---

**Nota:** Per dettagli sulle singole modifiche, consultare i file `nelder_mead_3A.c/h`, `nelder_mead_4A.c/h` (Crazyflie) e `nelder_mead_ORIG.c/h` (Matteo).
