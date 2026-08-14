"""Public three-dimensional axis-aligned bounding boxes."""
module BoundingBoxes

import StaticArrays: SVector

"""An axis-aligned bounding box with three-dimensional bounds."""
struct BoundingBox
    lower::SVector{3, Float64}
    upper::SVector{3, Float64}

    function BoundingBox(lower::SVector{3, Float64}, upper::SVector{3, Float64})
        all(lower .<= upper) || throw(ArgumentError("bounding-box lower bounds must not exceed upper bounds"))
        new(lower, upper)
    end
end

function BoundingBox(lower::AbstractVector, upper::AbstractVector)
    length(lower) == 3 || throw(DimensionMismatch("bounding boxes must be three-dimensional"))
    length(upper) == 3 || throw(DimensionMismatch("bounding boxes must be three-dimensional"))
    BoundingBox(SVector{3, Float64}(lower), SVector{3, Float64}(upper))
end

function BoundingBox(center::AbstractVector, radius::Number)
    length(center) == 3 || throw(DimensionMismatch("bounding boxes must be three-dimensional"))
    radius >= 0 || throw(ArgumentError("bounding-box radius must be nonnegative"))
    BoundingBox(center .- radius, center .+ radius)
end

function BoundingBox(radius::Number)
    radius >= 0 || throw(ArgumentError("bounding-box radius must be nonnegative"))
    BoundingBox(fill(-radius, 3), fill(radius, 3))
end

BoundingBox(box::BoundingBox) = box

lower(box::BoundingBox) = box.lower
upper(box::BoundingBox) = box.upper

end
