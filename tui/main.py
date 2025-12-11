#!/usr/bin/env python3

import sys
import os

# --- CONFIGURAZIONE PERCORSI SISTEMA ---
# 1. Ottiene il percorso della cartella corrente (tui/)
current_dir = os.path.dirname(os.path.abspath(__file__))
# 2. Ottiene il percorso della cartella genitore (PseudoBlast-Project/)
parent_dir = os.path.dirname(current_dir)
# 3. Aggiunge la root al path per poter importare 'src'
sys.path.append(parent_dir)

# --- IMPORT DALLA LOGICA CONDIVISA (SRC) ---
try:
    # Tenta di importare le funzioni dal modulo algoricerca dentro src
    from src.algoricerca import (
        sequenza_valida,
        choose_k,
        build_kmer_index,
        search_sequence,
        highlight_mismatches,
        save_csv,
        fuzzy_search,
        load_genomes,
        confronta_genomi
    )
except ImportError as e:
    print(f"\n[ERRORE CRITICO] Impossibile importare i moduli da 'src'.")
    print(f"Dettaglio errore: {e}")
    print(f"Verifica che esista: {os.path.join(parent_dir, 'src', 'algoricerca.py')}")
    sys.exit(1)

# --- IMPORT LIBRERIE ESTERNE ---
try:
    from Bio.Data import CodonTable
except ImportError:
    print("[ATTENZIONE] Biopython non installato. L'opzione 3 non funzionerà.")
    CodonTable = None

# --- CONFIGURAZIONE DATI ---
# Ora punta alla cartella centralizzata 'data/genomi'
ROOT = os.path.join(parent_dir, "data", "genomi")

def seleziona_specie_da_root():
    if not os.path.exists(ROOT):
        print(f"\n[ERRORE] La cartella dati non esiste: {ROOT}")
        print("Assicurati di aver creato la cartella 'data/genomi' nella root del progetto.")
        return None

    specie_disponibili = [
        d for d in os.listdir(ROOT) 
        if os.path.isdir(os.path.join(ROOT, d))
    ]
    specie_disponibili.sort()
    
    if not specie_disponibili:
        print(f" Nessuna sottocartella trovata dentro {ROOT}.")
        return None

    print(f"\n--- Specie disponibili in database ---")
    for i, nome_specie in enumerate(specie_disponibili, 1):
        print(f"{i}) {nome_specie}")
    
    while True:
        try:
            scelta_input = input("\nScegli il numero della specie: ").strip()
            idx = int(scelta_input)
            
            if 1 <= idx <= len(specie_disponibili):
                nome_scelto = specie_disponibili[idx - 1]
                path_completo = os.path.join(ROOT, nome_scelto)
                return path_completo
            else:
                print(f" Inserisci un numero tra 1 e {len(specie_disponibili)}.")
        except ValueError:
            print(" Input non valido. Inserisci un numero intero.")

def process_search_for_genome(genomes, query):
    """
    Esegue la ricerca della sequenza su un genoma caricato.
    """
    K = choose_k(len(query))
    print(f"\n K-mer automatico: K = {K}")
    print(" Costruzione indice...")
    index = build_kmer_index(genomes, K)

    max_mm = 0
    while True:
        try:
            inp = input("Numero massimo di mismatch (Invio per 0): ").strip()
            if inp == "":
                max_mm = 0
            else:
                max_mm = int(inp)
            break
        except ValueError:
            print("Inserisci un numero intero.")

    print("\n--- Ricerca esatta / mismatch ---\n")
    matches = search_sequence(query, index, genomes, K, max_mm)

    if matches:
        print(f" Trovati {len(matches)} risultati:\n")
        for m in matches:
            strand = "-" if m["rev"] else "+"
            print(f"- {m['file']} | {m['header']} | Pos: {m['pos']} | Strand: {strand} | Mismatches: {m['mismatches']}")
            print("  Segm:", highlight_mismatches(m["segment"], query))
            print("  Pos :", m["mm_positions"], "\n")

        if input("Vuoi salvare i risultati in CSV? (s/n): ").lower() == "s":
            save_csv(matches)
    else:
        print("Nessun match trovato.")

        if input("Vuoi provare una ricerca fuzzy (Smith-Waterman)? (s/n): ").lower() == "s":
            soglia = 20
            try:
                print("Lo score è la potenza di allineamento, più è alto, più le due sequenze sono simili.")
                inp = input("Soglia score (Invio per 20): ").strip()
                if inp: soglia = int(inp)
                
            except: pass
            
            print("\n Ricerca fuzzy a blocchi in corso...")
            
            fuzzy_results = fuzzy_search(query, genomes, soglia)
            
            if fuzzy_results:
                print(f"\nTrovati {len(fuzzy_results)} allineamenti fuzzy (Top 50):\n")
                
                for r in fuzzy_results[:50]:
                    # Recupera la posizione, usa 0 se non c'è
                    pos = r.get('pos_start', 0) 
                    
                    print(f"File: {r['file']} | Header: {r['header']}")
                    print(f"Posizione stimata: {pos} | Score: {r['score']}")
                    print(f"Query    : {r['q_aln']}")
                    print(f"Soggetto: {r['s_aln']}\n")
                # -------------------------
            else:
                print("Nessun risultato fuzzy rilevante.")

def main():
    print("\n" )
    print("############# PSEUDOBLAST (TUI) ##############")
    print(f"Database Genomi: {ROOT}")

    while True:
        print("\n" + "="*40)
        print("                  MENU")
        print("="*40)
        print("1) Cerca una sequenza in una specie")
        print("2) Confronta due specie (Genome vs Genome)")
        print("3) Stampa la tabella dei codoni")
        print("4) Esci")
        print("-" * 40)
        
        scelta = input("Opzione: ").strip()

        # -------------------------
        # OPZIONE 1: SEQUENZA -> SCELTA SPECIE
        # -------------------------
        if scelta == "1":
            query = input("\nInserisci la sequenza nucleotidica: ").strip().upper()
            if not sequenza_valida(query):
                print(" Sequenza non valida (usa solo A,C,G,T).")
                continue
            
            path_genoma = seleziona_specie_da_root()
            if not path_genoma: continue

            print(f"\n Caricamento genoma da: {os.path.basename(path_genoma)}...")
            genoma = load_genomes(path_genoma)
            
            if genoma:
                process_search_for_genome(genoma, query)
            else:
                print("Cartella vuota o nessun file FASTA trovato.")

        # -------------------------
        # OPZIONE 2: CONFRONTO DUE SPECIE
        # -------------------------
        elif scelta == "2":
            print("\nSELEZIONE PRIMA SPECIE:")
            path1 = seleziona_specie_da_root()
            if not path1: continue

            print("\nSELEZIONE SECONDA SPECIE:")
            path2 = seleziona_specie_da_root()
            if not path2: continue

            if path1 == path2:
                print(" Attenzione: hai selezionato la stessa specie due volte.")
                if input("Vuoi continuare comunque? (s/n): ").lower() != 's':
                    continue

            # Caricamento
            print("\n Caricamento Genoma 1...")
            genoma1 = load_genomes(path1)
            print(" Caricamento Genoma 2...")
            genoma2 = load_genomes(path2)

            if not genoma1 or not genoma2:
                print(" Uno dei due genomi non contiene sequenze valide.")
                continue

            soglia = 20
            try:
                inp = input("\nSoglia score minimo (Invio per 20): ").strip()
                if inp: soglia = int(inp)
            except: pass
        
            print(f"\n Confronto in corso (soglia {soglia})... attendere...")
            
            # Chiamata alla funzione (che ora usa i blocchi)
            compres = confronta_genomi(genoma1, genoma2, soglia)

            if compres:
                print(f" Trovati {len(compres)} allineamenti significativi (mostro i top 50):\n")
                
                for r in compres[:50]:
                    posA = r.get('posA_start', 0)
                    posB = r.get('posB_start', 0)
                    
                    print(f"Score: {r['score']} | {r['fileA']} (pos {posA}) <--> {r['fileB']} (pos {posB})")
                    print(f"   SeqA: {r['q_aln']}")
                    print(f"   SeqB: {r['s_aln']}\n")

            else:
                print("Nessun allineamento trovato.")
                
        # -------------------------
        # OPZIONE 3: STAMPO LA TABELLA DEI CODONI
        # -------------------------        
        elif scelta == "3": 
            if CodonTable:
                print("\n--- Tabella DNA Standard ---")
                print(CodonTable.standard_dna_table)
            else:
                print("Errore: Libreria BioPython non trovata.")
            continue
            
        # -------------------------
        # USCITA
        # -------------------------
        elif scelta == "4":
            print("Uscita.")
            break

        else:
            print("Scelta non valida.")

if __name__ == "__main__":
    main()
