function parse_section(section::PulseqSection{:adc}; shapes::Dict{Int, PulseqShape}, version::VersionNumber, kwargs...)
    result = Dict{Int, PulseqADC}()
    for line in section.content
        if version < v"1.5.0"
            props = parse_pulseq_dict(
                line,
                [:id, :num, :dwell, :delay, :freq, :phase],
                [Int, Int, Float64, Int, Float64, Float64],
            )
            result[props[:id]] = PulseqADC(
                props[:num],
                props[:dwell],
                props[:delay],
                props[:freq],
                props[:phase],
            ) 
        else
            props = parse_pulseq_dict(
                line,
                [:id, :num, :dwell, :delay, :freq_ppm, :phase_ppm, :freq, :phase, :phase_shape_id],
                [Int, Int, Float64, Int, Float64, Float64, Float64, Float64, Int],
            )
            result[props[:id]] = PulseqADC(
                props[:num],
                props[:dwell],
                props[:delay],
                props[:freq_ppm],
                props[:phase_ppm],
                props[:freq],
                props[:phase],
                _get_component(props[:phase_shape_id], shapes)
            ) 
        end
    end
    return result
end

function gen_section(comp:: PulseqComponents, ::Val{:adc})
    res = PulseqSection{:adc}(String[])
    for i in sort([keys(comp.adc)...])
        adc = comp.adc[i]
        values = string.(Any[
            i,
            adc.num,
            adc.dwell,
            adc.delay,
            adc.frequency_ppm,
            adc.phase_ppm,
            adc.frequency,
            adc.phase,
            add_components!(comp, adc.phase_shape)
        ])
        push!(res.content, join(values, " "))
    end
    return res
end