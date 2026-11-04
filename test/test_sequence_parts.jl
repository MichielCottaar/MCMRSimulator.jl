@testset "SequencePart promotion and MultSequencePart storage" begin
    parts = mr.SequenceParts
    zero_gradient = zero(SVector{3, Float64})
    gradient = SVector(1., 2., 3.)
    empty_part = parts.EmptyPart()
    constant_part = parts.ConstantPart(gradient)
    linear_part = parts.LinearPart(gradient, 2gradient)

    @test promote_type(parts.EmptyPart, parts.ConstantPart) == parts.ConstantPart
    @test promote_type(parts.EmptyPart, parts.LinearPart) == parts.LinearPart
    @test promote_type(parts.ConstantPart, parts.LinearPart) == parts.LinearPart
    @test convert(parts.ConstantPart, empty_part).strength == zero_gradient
    @test convert(parts.LinearPart, empty_part) == parts.LinearPart(zero_gradient, zero_gradient)
    @test convert(parts.LinearPart, constant_part) == parts.LinearPart(gradient, gradient)

    pulse = parts.ConstantPulse(0.25, 30., 1.)
    pulse_empty = parts.PulsePart([pulse], empty_part, 0.4, 0.8)
    pulse_linear = parts.PulsePart([pulse], linear_part, 0.5, 0.7)
    @test promote_type(typeof(pulse_empty), typeof(pulse_linear)) == typeof(pulse_linear)
    promoted_pulse = convert(typeof(pulse_linear), pulse_empty)
    @test promoted_pulse.gradient == convert(parts.LinearPart, empty_part)
    @test promoted_pulse.pulse === pulse_empty.pulse
    @test promoted_pulse.first_duration == pulse_empty.first_duration
    @test promoted_pulse.last_duration == pulse_empty.last_duration

    function mult_part(seq_parts; VT=SVector{length(seq_parts)})
        return parts.MultSequencePart{VT}(1., seq_parts, fill(nothing, length(seq_parts)))
    end

    no_rf = mult_part([empty_part, constant_part, linear_part])
    @test eltype(no_rf.parts) == parts.LinearPart
    @test no_rf.parts == [convert(parts.LinearPart, empty_part), convert(parts.LinearPart, constant_part), linear_part]

    rf = mult_part([pulse_empty, pulse_linear])
    @test eltype(rf.parts) == typeof(pulse_linear)
    @test rf.parts[1].gradient == convert(parts.LinearPart, empty_part)

    mixed = mult_part([empty_part, parts.PulsePart([pulse], linear_part)])
    @test eltype(mixed.parts) == Union{parts.EmptyPart, typeof(parts.PulsePart([pulse], linear_part))}
    @test mixed.parts[1] == empty_part
    @test mixed.parts[2] isa parts.PulsePart{parts.LinearPart}

    large_mixed = mult_part(vcat(fill(empty_part, 25), fill(pulse_empty, 25)); VT=Vector)
    @test eltype(large_mixed.parts) == Union{parts.EmptyPart, typeof(pulse_empty)}
    @test length(large_mixed.parts) == 50
    @test all(p -> p isa parts.EmptyPart, large_mixed.parts[1:25])
end
