import os
import sys

# --- 1. CONFIGURAZIONE PERCORSI (PRIMA DI TUTTO IL RESTO) ---
# Ottieni la cartella corrente (web/)
current_dir = os.path.dirname(os.path.abspath(__file__))
# Ottieni la cartella genitore (PseudoBlast-Project/)
parent_dir = os.path.dirname(current_dir)
# Aggiungi la root al percorso di sistema
sys.path.append(parent_dir)

# --- 2. IMPORT ---
from flask import Flask, render_template, request

# Ora possiamo importare 'src' perché sys.path è aggiornato
try:
    from src.algoricerca import load_genomes, confronta_genomi, fuzzy_search
except ImportError as e:
    print(f"ERRORE CRITICO: Non trovo il modulo 'src'. {e}")
    sys.exit(1)

app = Flask(__name__)

# --- 3. CONFIGURAZIONE DATI ---
# Aggiornato per puntare alla cartella centralizzata 'data/genomi'
ROOT = os.path.join(parent_dir, "data", "genomi")

def get_species_list():
    """Trova le cartelle delle specie disponibili."""
    if not os.path.exists(ROOT): 
        print(f"Attenzione: La cartella dati non esiste: {ROOT}")
        return []
    
    # Filtra solo le cartelle
    species = [
        d for d in os.listdir(ROOT) 
        if os.path.isdir(os.path.join(ROOT, d))
    ]
    return sorted(species)

@app.route('/', methods=['GET', 'POST'])
def home():
    specie = get_species_list()
    results = []
    message = ""
    
    if request.method == 'POST':
        action = request.form.get('action')

        # Parametri opzionali
        try:
            min_score = int(request.form.get("min_score", 0))
        except ValueError:
            min_score = 0

        # --- CASO 1: RICERCA SEQUENZA ---
        if action == 'search_seq':
            seq = request.form.get('sequence', '').upper().strip()
            specie_scelta = request.form.get('genome_single')
            
            if seq and specie_scelta:
                path = os.path.join(ROOT, specie_scelta)
                genoma = load_genomes(path)

                if genoma:
                    # CORREZIONE QUI: cambiato 'score_threshold' in 'soglia'
                    try:
                        raw_results = fuzzy_search(
                            seq,
                            genoma,
                            soglia=min_score 
                        )
                        # Se ti dà ancora errore, prova a togliere "soglia=" e passare solo min_score
                    except TypeError:
                         # Fallback: prova a passarlo come argomento posizionale se il nome è diverso
                         raw_results = fuzzy_search(seq, genoma, min_score)

                    results = raw_results[:50] if raw_results else []
                    message = f"Ricerca completata su {specie_scelta}: {len(results)} risultati trovati."
                else:
                    message = "Errore nel caricamento del genoma."

        # --- CASO 2: CONFRONTO GENOMI ---
        elif action == 'compare_genomes':
            sp1 = request.form.get('genome_1')
            sp2 = request.form.get('genome_2')
            
            if sp1 and sp2:
                gen1 = load_genomes(os.path.join(ROOT, sp1))
                gen2 = load_genomes(os.path.join(ROOT, sp2))
                
                if gen1 and gen2:
                    # CORREZIONE QUI: cambiato 'score_threshold' in 'soglia'
                    try:
                        raw_results = confronta_genomi(
                            gen1,
                            gen2,
                            soglia=min_score
                        )
                    except TypeError:
                        # Fallback posizionale
                        raw_results = confronta_genomi(gen1, gen2, min_score)

                    results = raw_results[:50] if raw_results else []
                    message = f"Confronto tra {sp1} e {sp2} completato: {len(results)} match trovati."
                else:
                    message = "Errore nel caricamento dei genomi."

    return render_template('index.html', specie=specie, results=results, message=message)

if __name__ == '__main__':
    # host='0.0.0.0' permette di vedere il sito dagli altri computer nella rete locale
    app.run(debug=True, host='0.0.0.0', port=5000)
