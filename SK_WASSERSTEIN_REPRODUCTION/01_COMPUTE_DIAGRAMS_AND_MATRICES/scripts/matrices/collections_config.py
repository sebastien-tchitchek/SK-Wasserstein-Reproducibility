#!/usr/bin/env python3

from __future__ import annotations

COLLECTIONS = (
    {"key": "2004_isabel_3D", "label": "Isabel 3D", "paper_label": "Isabel", "expected": 12, "k": 3},
    {"key": "2006_earthquake_3D", "label": "Earthquake 3D", "paper_label": "Earthquake", "expected": 12, "k": 3},
    {"key": "2008_ionization_front_2D", "label": "Ionization Front 2D", "paper_label": "Ionization 2D", "expected": 16, "k": 4},
    {"key": "2008_ionization_front_3D", "label": "Ionization Front 3D", "paper_label": "Ionization 3D", "expected": 16, "k": 4},
    {"key": "2014_volcanic_eruptions_2D", "label": "Volcanic Eruptions 2D", "paper_label": "Volcanic", "expected": 12, "k": 3},
    {"key": "2016_viscous_fingering_3D", "label": "Viscous Fingering 3D", "paper_label": "Viscous Fingering", "expected": 15, "k": 3},
    {"key": "2017_cloud_processes_2D", "label": "Cloud Processes 2D", "paper_label": "Cloud Processes", "expected": 12, "k": 3},
    {"key": "2018_asteroid_impact_3D_clustering", "label": "Asteroid Impact 3D Clustering", "paper_label": "Asteroid ensemble", "expected": 7, "k": 2},
    {"key": "2018_asteroid_impact_3D_temporal_subsampling", "label": "Asteroid Impact 3D Temporal Subsampling", "paper_label": "Asteroid temporal", "expected": 20, "k": 4},
    {"key": "sea_surface_height", "label": "Sea Surface Height 2D", "paper_label": "Sea Surface Height", "expected": 48, "k": 4},
    {"key": "starting_vortex", "label": "Starting Vortex 2D", "paper_label": "Starting Vortex", "expected": 12, "k": 2},
    {"key": "vortex_street", "label": "Vortex Street 2D", "paper_label": "Vortex Street", "expected": 45, "k": 5},
)

DEFAULT_LEVELS = (10, 20, 30, 40)
PRIMARY_LEVEL = 30
REFERENCE_LEVEL = 40
EXPECTED_TOTAL_DIAGRAMS = sum(int(c["expected"]) for c in COLLECTIONS)
EXPECTED_MATRICES_PER_COLLECTION = 2 * len(DEFAULT_LEVELS) + 1
EXPECTED_TOTAL_MATRICES = len(COLLECTIONS) * EXPECTED_MATRICES_PER_COLLECTION
EXPECTED_PAIRWISE_COMPARISONS = sum(
    int(c["expected"]) * (int(c["expected"]) - 1) // 2 for c in COLLECTIONS
)
EXPECTED_PERSISTENCE_PAIRS = 3_376_524
