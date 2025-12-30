import math

def calculate_binding_energy(motif_count, g_motif=-8.2, coop_factor=0.15):
    """
    Calculates the cumulative Gibbs Free Energy (Delta G) for a cluster 
    of transcription factor binding sites.
    
    Parameters:
    motif_count (int): Number of tandem repeats (e.g., 8 for Bowhead).
    g_motif (float): Individual binding energy of the motif in kcal/mol.
    coop_factor (float): Cooperative binding energy per adjacent site.
    
    Returns:
    float: Total Delta G in kcal/mol.
    """
    if motif_count == 0:
        return 0
    
    # Base energy for the motifs
    total_g = motif_count * g_motif
    
    # Energy contribution from cooperative interactions between adjacent sites
    # Total adjacent pairs = motif_count - 1
    total_coop = (motif_count - 1) * coop_factor if motif_count > 1 else 0
    
    # In this model, high recruitment (super-enhancement) is 
    # represented by the additive effect of the sites.
    # Note: We use the cumulative sum for the 'volume knob' effect.
    return round(total_g + total_coop, 2)

# Test Data based on Discovery Locus (Bmy-CIRBP-TR8)
species_data = {
    "Human": 1,
    "Porpoise": 1,
    "Minke": 2,
    "Bowhead": 8,
    "Greenland Shark": 14
}

print("Species | Repeats | Predicted Delta G (kcal/mol)")
print("-" * 45)
for species, count in species_data.items():
    delta_g = calculate_binding_energy(count)
    print(f"{species:15} | {count:7} | {delta_g:12}")
