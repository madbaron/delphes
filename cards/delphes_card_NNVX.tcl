#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  PileUpMerger
  ParticlePropagator

  TrackMergerAll
  SiliconHits
  DisplacedVertexDumper

  TreeWriter
}

###############
# PileUp Merger
###############

module PileUpMerger PileUpMerger {
  set InputArray Delphes/stableParticles

  set ParticleOutputArray stableParticles
  set VertexOutputArray vertices

  # pre-generated minbias input file
  set PileUpFile /Applications/Delphes/MinBias.pileup

  # average expected pile up
  # set MeanPileUp 60
  set MeanPileUp 0

  # maximum spread in the beam direction in m
  set ZVertexSpread 0.075
  
  # maximum spread in time in s
  set TVertexSpread 800E-12

  # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)
  #set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}
  set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))}

}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray PileUpMerger/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.15
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.51

  # magnetic field
  set Bz 2.0
}

##############
# Track merger all
##############

module Merger TrackMergerAll {
# add InputArray InputArray
  add InputArray ParticlePropagator/chargedHadrons
  add InputArray ParticlePropagator/electrons
  add InputArray ParticlePropagator/muons
  set OutputArray tracks
}

##############
# Displaced Vertices
##############

module DisplacedVertexDumper DisplacedVertexDumper {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  set InputArray TrackMergerAll/tracks

}

##############
# Silicon Hits
##############

module SiliconHits SiliconHits {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  set InputArray TrackMergerAll/tracks

}

##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {

  add Branch TrackMergerAll/tracks Track Track
  add Branch PileUpMerger/vertices Vertex Vertex
  add Branch SiliconHits/hits Hits Hit 
  add Branch DisplacedVertexDumper/displacedvertices DisplacedVertices DisplacedVertex 

}
