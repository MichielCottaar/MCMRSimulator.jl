"""
    ObstructionProperties(kwargs...)
    ObstructionProperties(obstruction)

The generic properties of an obstruction that all obstructions have in common. They include:
- `unique_identifier`: unique identifier
- `MT_fraction`: fraction of signal transfered in case of a collision
"""
struct ObstructionProperties
    id :: UUID
    MT_fraction :: Float
    permeability :: Float
end
ObstructionProperties(; MT_fraction=zero(Float), permeability=zero(Float)) = ObstructionProperties(uuid1(), MT_fraction, permeability)

for accessor in [:id, :MT_fraction]
    @eval $(accessor)(prop :: ObstructionProperties) = prop.$(accessor)
    @eval $(accessor)(obstruction) = $(accessor)(ObstructionProperties(obstruction))
end