module InstantReadouts
import ...BuildingBlocks: BuildingBlock, BuildingBlockPlaceholder, duration, properties, to_mcmr_components
import ...SequenceBuilders: SequenceBuilder, start_time

struct InstantReadout <: BuildingBlock
    builder::SequenceBuilder
end

InstantReadout() = BuildingBlockPlaceholder{InstantReadout}()

duration(::InstantReadout) = 0.
properties(::Type{<:InstantReadout}) = []

to_mcmr_components(readout::InstantReadout) = Readout(value(start_time(readout)))
end