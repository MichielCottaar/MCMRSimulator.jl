using CairoMakie
using FFTW
using LinearAlgebra

sz = 2^7
half = Int(sz/2)

function fill_chi(sz)
    chi = zeros(sz, sz, sz);
    for index in CartesianIndices(chi)
        radius = norm(Tuple(index) .- [half, half, half])
        if radius <= 30
            chi[index] = 1.
        end
    end
    chi
end
chi = fill_chi(sz)

image(chi[:, :, half]);

transformed = fft(chi);
image(abs.(transformed[:, :, half]))


function get_mult(sz)
    k_single = fftfreq(sz)
    kx = reshape(k_single, sz, 1, 1)
    ky = reshape(k_single, 1, sz, 1)
    kz = reshape(k_single, 1, 1, sz)
    @. 1/3 - kz * kz / (kx * kx + ky * ky + kz * kz + 1e-12)
end
to_mult = get_mult(sz);


field = ifft(to_mult .* transformed);
image(to_mult[:, half, :])

image(real.(field[:, half, :]))

lines(real.(field[:, half, half]))
lines(real.(field[half, half, :]))