"""
    Annulus(inner_radius, outer_radius; chi_I=-0.1, chi_A=-0.1, myelin=false)

Create an annulus with an inner and outer radius.
Using [`annuli`](@ref) is recommended.
Water can freely diffuse within the inner radius, and between the inner and outer radius, but not in between.
Myelin can be added by setting `myelin` to true.
Restrictions are only present at the inner and outer cylinder, not in the myelin in between.
"""
struct Annulus <: BaseObstruction{2}
    inner :: Cylinder
    outer :: Cylinder
    myelin :: Bool
    chi_I :: Float
    chi_A :: Float
    internal_field :: Float
    external_field :: Float
    function Annulus(inner::Real, outer::Real; chi_I=-0.1, chi_A=-0.1, myelin=false)
        @assert outer > inner
        inner = Float(inner)
        outer = Float(outer)
        internal_field = 0.75 * chi_A * log(outer / inner)
        external_field = (chi_I + chi_A / 4) * (outer * outer - inner * inner) / 2
        new(Cylinder(inner), Cylinder(outer), myelin, chi_I, chi_A, internal_field, external_field)
    end
end

Base.copy(a::Annulus) = Annulus(a.inner.radius, a.outer.radius; chi_I=a.chi_I, chi_A=a.chi_A, myelin=a.myelin)
isinside(a::Annulus, pos::SVector{2, Float}) = isinside(a.inner, pos) + isinside(a.outer, pos)
BoundingBox(a::Annulus) = BoundingBox(a.outer)

function total_susceptibility(a::Annulus)
    if !a.myelin
        return zero(Float)
    end
    chi = a.chi_I + a.chi_A / 4
    2 * π * chi * (a.outer.radius^2 - a.inner.radius^2)
end

"""
    annuli(inner, outer; myelin=false, chi_I=-0.1, chi_A=-0.1, positions=[0, 0], repeats=[Inf, Inf], rotation=I(3)

Creates one or more [`Annulus`](@ref) with given `inner` and `outer` radii.
[Myelinated annuli](@ref) can be created by setting the `myelin` to true.
All parameters can be either a single value or a vector of values.

The `positions`, `repeats`, and `rotation` control the annulus position and orientation and is explained in 
more detail in [Defining the goemetry](@ref).
"""
function annuli(args...; kwargs...)
    TransformObstruction(Annulus, args...; kwargs...)
end

function detect_collision(movement :: Movement{2}, annulus :: Annulus, previous=empty_collision)
    if previous !== empty_collision
        if previous.id == annulus.inner.id
            if previous.index == 1
                # we are inside the inner sphere
                return detect_collision(movement, annulus.inner, previous)
            else
                # we just bounced on the outside of the inner sphere. Let's not do that again
                return detect_collision(movement, annulus.outer, previous)
            end
        elseif previous.id == annulus.outer.id
            if previous.index == 0
                # We just bounced on the outside of the outer cylinder. Let's check if we hit it again (could be due to repeats).
                return detect_collision(movement, annulus.outer, previous)
            end
        end
    end
    # Not sure where we are...
    c_inner = detect_collision(movement, annulus.inner, previous)
    c_outer = detect_collision(movement, annulus.outer, previous)
    c_inner.distance < c_outer.distance ? c_inner : c_outer
end


function lorentz_off_resonance(annulus::Annulus, position::SVector{2, Float}, b0_field::SVector{2, Float}, repeat_dist::SVector{2, Float}, radius::Float, nrepeats::SVector{2, Int})
    field = zero(Float)
    if !annulus.myelin
        return field
    end
    lorentz_radius_sq = radius * radius
    sin_theta_sq = b0_field[1] * b0_field[1] - b0_field[2] + b0_field[2]
    if iszero(sin_theta_sq)
        return field
    end
    for i in -nrepeats[1]:nrepeats[1]
        xshift = i * repeat_dist[1]
        p1 = position[1] + xshift
        for j in -nrepeats[2]:nrepeats[2]
            yshift = j * repeat_dist[2]
            p2 = position[2] + yshift
            rsq = p1 * p1 + p2 * p2
            if rsq > lorentz_radius_sq
                continue
            end
            if i == 0 && j == 0 && rsq < annulus.inner.radius^2
                # inside inner cylinder
                field += annulus.internal_field * sin_theta_sq
            elseif i == 0 && j == 0 && rsq < annulus.outer.radius^2
                # between cylinders
                cos2 = (b0_field[1] * p1 + b0_field[2] * p2)^2 / rsq
                cos2f = 2 * cos2 - 1
                field += annulus.chi_I * (2//3 -  sin_theta_sq * (1 + cos2f * (annulus.inner.radius^2 / rsq))) / 2
                field += annulus.chi_A * sin_theta_sq * (-5//12 - cos2f/8 * (1 + annulus.inner.radius^2/rsq) + 3//8 * log(annulus.outer.radius^2/rsq)) - (1-sin_theta_sq) / 6
            else
                # outside annulus
                cos2 = (b0_field[1] * p1 + b0_field[2] * p2)^2 / rsq
                field += annulus.external_field * sin_theta_sq * (2 * cos2 - 1) / rsq
            end
        end
    end
    return field
end