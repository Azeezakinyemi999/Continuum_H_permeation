"""Hydrogen permeation model, Levels 1-6.

A hierarchy of steady-state analytical models for hydrogen transport through an
oxide-coated alloy wall: perfect metal (L1), oxide layer (L2), defective oxide
(L3), defective metal with grain-boundary enhancement and trapping (L4), the full
defective system (L5), and surface dissociation kinetics (L6).

All material data and operating conditions resolve through the single switch in
:mod:`calculations.config.model_config`; read that module before changing a study.

This package deliberately re-exports nothing. Import from the module that owns the
function, so each name has exactly one import path.
"""
