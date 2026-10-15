# ==============================================================================
# Analysis 2: Boolean Network Dynamics Modeling (BoolNet)
# Translates the dual/three-phase mechanistic hypothesis into a synchronous
# Boolean network and tests: (1) signal propagation from PBM input,
# (2) IL6 vs ACTB knockout dynamics, (3) hysteresis/self-sustaining behavior
# ==============================================================================

suppressMessages(library(BoolNet))

net <- loadNetwork("scripts/pbm_boolean_network.txt")
cat("=== Network loaded ===\n")
print(net)
cat("\n")

all_genes <- net$genes
cat("Nodes:", paste(all_genes, collapse=", "), "\n\n")

# ------------------------------------------------------------------------------
# Helper: run a synchronous trajectory from a given start state for k steps
# ------------------------------------------------------------------------------
run_trajectory <- function(network, start_state, steps) {
  traj <- matrix(NA, nrow = steps + 1, ncol = length(network$genes))
  colnames(traj) <- network$genes
  traj[1, ] <- start_state
  state <- start_state
  for (i in 1:steps) {
    state <- stateTransition(network, state)
    traj[i + 1, ] <- state
  }
  traj
}

# ==============================================================================
# SIMULATION 1: Baseline -- PBM OFF (negative control)
# ==============================================================================
cat("=================================================================\n")
cat("SIMULATION 1: PBM = OFF (negative control)\n")
cat("=================================================================\n")
net_off <- fixGenes(net, "PBM", 0)
start0 <- setNames(rep(0, length(all_genes)), all_genes)
traj_off <- run_trajectory(net_off, start0, 8)
print(traj_off)
cat("Final state all zero:", all(traj_off[nrow(traj_off), ] == 0), "\n\n")

# ==============================================================================
# SIMULATION 2: PBM = ON, full network (test propagation to ACTB/CTNNB1/Drainage)
# ==============================================================================
cat("=================================================================\n")
cat("SIMULATION 2: PBM = ON (full network, wild-type)\n")
cat("=================================================================\n")
net_on <- fixGenes(net, "PBM", 1)
start1 <- setNames(rep(0, length(all_genes)), all_genes)
start1["PBM"] <- 1
traj_on <- run_trajectory(net_on, start1, 14)
print(traj_on)
write.csv(traj_on, "results/boolean_trajectory_wildtype.csv", row.names = FALSE)

final_wt <- traj_on[nrow(traj_on), ]
cat("\nFinal (attractor) state:\n")
print(final_wt)
cat("\nACTB activated:", final_wt["ACTB"] == 1, "\n")
cat("CTNNB1 activated:", final_wt["CTNNB1"] == 1, "\n")
cat("Drainage activated:", final_wt["Drainage"] == 1, "\n\n")

# Formal attractor detection
attr_wt <- getAttractors(net_on)
cat("Number of attractors (PBM fixed ON):", length(attr_wt$attractors), "\n")
print(getAttractorSequence(attr_wt, 1))
cat("\n")

# ==============================================================================
# SIMULATION 3: IL6 knockout (PBM=ON, IL6 fixed OFF) -- dynamic validation
# Run to full attractor length (14 steps, chain depth is ~11)
# ==============================================================================
cat("=================================================================\n")
cat("SIMULATION 3: IL6 knockout (PBM=ON, IL6=OFF)\n")
cat("=================================================================\n")
net_il6ko <- fixGenes(net, c("PBM", "IL6"), c(1, 0))
traj_il6ko <- run_trajectory(net_il6ko, start1, 14)
print(traj_il6ko)
write.csv(traj_il6ko, "results/boolean_trajectory_IL6ko.csv", row.names = FALSE)
final_il6ko <- traj_il6ko[nrow(traj_il6ko), ]
cat("\nFinal state with IL6 knocked out:\n")
print(final_il6ko)
cat("\nACTB still activated (via TNF redundancy):", final_il6ko["ACTB"] == 1, "\n")
cat("Drainage still activated despite IL6 knockout:", final_il6ko["Drainage"] == 1, "\n")

attr_il6ko <- getAttractors(net_il6ko)
cat("Number of attractors (IL6 knockout):", length(attr_il6ko$attractors), "\n")
print(getAttractorSequence(attr_il6ko, 1))
cat("\n")

# ==============================================================================
# SIMULATION 4: ACTB knockout (PBM=ON, ACTB fixed OFF) -- structural bottleneck test
# ==============================================================================
cat("=================================================================\n")
cat("SIMULATION 4: ACTB knockout (PBM=ON, ACTB=OFF)\n")
cat("=================================================================\n")
net_actbko <- fixGenes(net, c("PBM", "ACTB"), c(1, 0))
traj_actbko <- run_trajectory(net_actbko, start1, 14)
print(traj_actbko)
write.csv(traj_actbko, "results/boolean_trajectory_ACTBko.csv", row.names = FALSE)
final_actbko <- traj_actbko[nrow(traj_actbko), ]
cat("\nFinal state with ACTB knocked out:\n")
print(final_actbko)
cat("\nCTNNB1 activated:", final_actbko["CTNNB1"] == 1, "\n")
cat("Drainage activated:", final_actbko["Drainage"] == 1, "\n")
cat("(Expected: FULL blockade, since CTNNB1 has no redundant input in this model)\n\n")

# ==============================================================================
# SIMULATION 5: Hysteresis test -- PBM pulses ON for 3 steps then OFF permanently.
# Does the KLF2 feedback loop sustain the activated state (self-sustaining Phase 3)?
# ==============================================================================
cat("=================================================================\n")
cat("SIMULATION 5: Transient PBM pulse (ON for 3 steps, then OFF) -- hysteresis test\n")
cat("=================================================================\n")
# Manual mixed simulation: PBM is NOT fixed here; we directly overwrite its value
# for the first 3 steps, then let the network evolve with PBM free (governed by
# its own rule PBM=PBM, i.e. it holds whatever value we last set unless overwritten)
state <- setNames(rep(0, length(all_genes)), all_genes)
state["PBM"] <- 1
pulse_steps <- 3
total_steps <- 14
traj_pulse <- matrix(NA, nrow = total_steps + 1, ncol = length(all_genes))
colnames(traj_pulse) <- all_genes
traj_pulse[1, ] <- state

for (i in 1:total_steps) {
  state <- stateTransition(net, state)
  if (i < pulse_steps) {
    state["PBM"] <- 1   # hold PBM on during the pulse window
  } else {
    state["PBM"] <- 0   # force PBM off from step `pulse_steps` onward
  }
  traj_pulse[i + 1, ] <- state
}
print(traj_pulse)
write.csv(traj_pulse, "results/boolean_trajectory_pulse_hysteresis.csv", row.names = FALSE)

final_pulse <- traj_pulse[nrow(traj_pulse), ]
cat("\nFinal state after PBM pulse removed (steady state):\n")
print(final_pulse)
cat("\nDrainage SUSTAINED after PBM removed:", final_pulse["Drainage"] == 1, "\n")
cat("(Tests whether the KLF2-\u003eACTB feedback loop creates a bistable/self-sustaining state)\n\n")

# Formal check: is the pulse trajectory converging to a fixed point or cycling?
tail_states <- traj_pulse[(nrow(traj_pulse)-4):nrow(traj_pulse), c("ACTB","CTNNB1","Permeability","ShearStress","KLF2","Drainage")]
cat("Last 5 states (ACTB / CTNNB1 / Permeability / ShearStress / KLF2 / Drainage):\n")
print(tail_states)
is_cycling <- !all(apply(tail_states, 2, function(x) length(unique(x))==1))
cat("\nSystem is CYCLING, not converged to a fixed point:", is_cycling, "\n\n")

# ==============================================================================
# SIMULATION 6: Memory-augmented model
# Adds a biologically-motivated persistence term to ACTB, representing the
# structural persistence of remodeled actin filaments (slow depolymerization
# kinetics) rather than instantaneous reversibility. Tests whether this
# resolves the oscillation into genuine bistable self-sustainment, as
# narratively proposed for Phase 3 of the manuscript's mechanistic model.
# ==============================================================================
cat("=================================================================\n")
cat("SIMULATION 6: Memory-augmented model (ACTB persistence term added)\n")
cat("=================================================================\n")

net_lines <- readLines("scripts/pbm_boolean_network.txt")
net_lines_mem <- gsub("^ACTB, \\(STAT3 \\| TNF\\) \\| KLF2$",
                       "ACTB, (STAT3 | TNF) | KLF2 | ACTB", net_lines)
writeLines(net_lines_mem, "scripts/pbm_boolean_network_memory.txt")
net_mem <- loadNetwork("scripts/pbm_boolean_network_memory.txt")
cat("Modified rule for ACTB now includes self-persistence term (see network file)\n\n")

state <- setNames(rep(0, length(all_genes)), all_genes)
state["PBM"] <- 1
traj_mem <- matrix(NA, nrow = total_steps + 1, ncol = length(all_genes))
colnames(traj_mem) <- all_genes
traj_mem[1, ] <- state
for (i in 1:total_steps) {
  state <- stateTransition(net_mem, state)
  if (i < pulse_steps) {
    state["PBM"] <- 1
  } else {
    state["PBM"] <- 0
  }
  traj_mem[i + 1, ] <- state
}
print(traj_mem)
write.csv(traj_mem, "results/boolean_trajectory_pulse_memory.csv", row.names = FALSE)

tail_mem <- traj_mem[(nrow(traj_mem)-4):nrow(traj_mem), c("ACTB","CTNNB1","Permeability","ShearStress","KLF2","Drainage")]
cat("Last 5 states, memory-augmented model:\n")
print(tail_mem)
is_stable_mem <- all(apply(tail_mem, 2, function(x) length(unique(x))==1))
cat("\nMemory-augmented model reaches a STABLE fixed point:", is_stable_mem, "\n")
final_mem <- traj_mem[nrow(traj_mem), ]
cat("Drainage sustained (memory model, PBM long removed):", final_mem["Drainage"] == 1, "\n\n")

cat("=================================================================\n")
cat("ALL BOOLEAN SIMULATIONS COMPLETE\n")
cat("=================================================================\n")
