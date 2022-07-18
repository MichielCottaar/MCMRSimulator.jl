function norm_angle(angle)
    angle = mod(angle, 360)
    if angle > 180
        angle -= 360
    end
    angle
end

struct SpinOrientation
    longitudinal :: Real
    transverse :: Real
    phase :: Real
end

struct Spin
    time :: Real
    position :: SVector{3,Real}
    orientation :: SpinOrientation
end
Spin(;time=0., position=zero(SVector{3,Real}), longitudinal=1., transverse=0., phase=0.) = Spin(time, position, SpinOrientation(longitudinal, transverse, deg2rad(phase)))

for param in (:longitudinal, :transverse)
    @eval $param(o :: SpinOrientation) = o.$param
    @eval $param(s :: Spin) = $param(s.orientation)
end
phase(o :: SpinOrientation) = norm_angle(rad2deg(o.phase))

phase(s :: Spin) = phase(s.orientation)
time(s :: Spin) = s.time
position(s :: Spin) = s.position