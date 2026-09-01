#pragma once

#include <vector>
#include "PersistenceDiagramUtils.h"

void diagram_to_v_tmid(const ttk::DiagramType &D,
                       int L,
                       std::vector<double> &out_v);

std::vector<std::vector<double>>
precompute_all_v(const std::vector<ttk::DiagramType> &geodesic, int L);

double pot_cost_on_SK_1d_from_sorted_tmid(const std::vector<double> &v1_sorted,
                                          const std::vector<double> &v2_sorted,
                                          int L,
                                          int choiceHilbertDistance, int GLevel);

void diagram_to_skot_sorted_vectors(const ttk::DiagramType &D,
                                    int L,
                                    std::vector<double> &out_off_sorted,
                                    std::vector<double> &out_diag_sorted);

double skot_distance_from_precomputed_sorted(
  const std::vector<double> &off1_sorted,
  const std::vector<double> &diag1_sorted,
  const std::vector<double> &off2_sorted,
  const std::vector<double> &diag2_sorted);

double skot_distance_normalized(const ttk::DiagramType &D1,
                                const ttk::DiagramType &D2,
                                int L);

double dsk_distance_normalized(const ttk::DiagramType &D1,
                               const ttk::DiagramType &D2,
                               int L);

void precompute_all_skot_sorted_vectors(
  const std::vector<ttk::DiagramType> &diagrams,
  int L,
  std::vector<std::vector<double>> &all_off_sorted,
  std::vector<std::vector<double>> &all_diag_sorted);

struct SkOffAtom {
  double t;
  double x;
  double y;
};

void diagram_to_w2_delta_sk_sorted_data(
  const ttk::DiagramType &D,
  int L,
  std::vector<SkOffAtom> &out_off_sorted,
  std::vector<double> &out_diag_sorted);

double w2_delta_sk_squared_from_precomputed_sorted(
  const std::vector<SkOffAtom> &off1_sorted,
  const std::vector<double> &diag1_sorted,
  const std::vector<SkOffAtom> &off2_sorted,
  const std::vector<double> &diag2_sorted);

double w2_delta_sk_from_precomputed_sorted(
  const std::vector<SkOffAtom> &off1_sorted,
  const std::vector<double> &diag1_sorted,
  const std::vector<SkOffAtom> &off2_sorted,
  const std::vector<double> &diag2_sorted);

double w2_delta_sk_squared_normalized(
  const ttk::DiagramType &D1,
  const ttk::DiagramType &D2,
  int L);

double w2_delta_sk_normalized(
  const ttk::DiagramType &D1,
  const ttk::DiagramType &D2,
  int L);

void precompute_all_w2_delta_sk_sorted_data(
  const std::vector<ttk::DiagramType> &diagrams,
  int L,
  std::vector<std::vector<SkOffAtom>> &all_off_sorted,
  std::vector<std::vector<double>> &all_diag_sorted);

double skot_distance_to_empty_from_precomputed_sorted(
  const std::vector<double> &off_sorted,
  const std::vector<double> &diag_sorted);

double dsk_distance_to_empty_from_precomputed_sorted(
  const std::vector<double> &off_sorted,
  const std::vector<double> &diag_sorted);

double w2_delta_sk_squared_to_empty_from_precomputed_sorted(
  const std::vector<SkOffAtom> &off_sorted,
  const std::vector<double> &diag_sorted);

double w2_delta_sk_to_empty_from_precomputed_sorted(
  const std::vector<SkOffAtom> &off_sorted,
  const std::vector<double> &diag_sorted);

double skot_distance_to_empty_normalized(
  const ttk::DiagramType &D,
  int L);

double dsk_distance_to_empty_normalized(
  const ttk::DiagramType &D,
  int L);

double w2_delta_sk_squared_to_empty_normalized(
  const ttk::DiagramType &D,
  int L);

double w2_delta_sk_to_empty_normalized(
  const ttk::DiagramType &D,
  int L);

void precompute_all_skot_distances_from_empty(const std::vector<ttk::DiagramType> &diagrams,
                                              int L,
                                              std::vector<double> &distances, 
                                              std::vector<std::vector<double>> &all_off_sorted,
                                              std::vector<std::vector<double>> &all_diag_sorted);

void precompute_all_w2_delta_sk_from_empty(const std::vector<ttk::DiagramType> &diagrams,
                                              int L,
                                              std::vector<double> &distances, 
                                              std::vector<std::vector<SkOffAtom>> &all_off_sorted,
                                              std::vector<std::vector<double>> &all_diag_sorted);


void precompute_all_pot_distances_from_empty(const std::vector<std::vector<double>> &allV,
                                              int L,
                                              int choiceHilbertDistance,
                                              std::vector<double> &distances, int GLevel);