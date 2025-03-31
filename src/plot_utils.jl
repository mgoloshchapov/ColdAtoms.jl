#Function for visualisation of samples
function samples_visualise_rus(samples; samples_harmonic=[], harmonic=false)
    x, y, z, vx, vy, vz = invert(samples);
    if harmonic
        xh, yh, zh, vxh, vyh, vzh = invert(samples_harmonic);
    end;

    figure(figsize=(8,8))
    subplot(221)
    hist2D(x, z, bins=50, range=[[-1.0, 1.0], [-2.5, 2.5]],cmap="plasma", rasterized=true);
    xlabel("x, мкм")
    ylabel("z, мкм")
    title("Распределение по координатам");

    subplot(222)
    hist(x, bins=[-1.5:0.02:1.5;], density=true, alpha=0.5, label="x")
    hist(z, bins=[-2.0:0.1:2.0;], density=true, alpha=0.5, label="z")
    if harmonic
        hist(xh, bins=[-1.5:0.01:1.5;], density=true, alpha=0.5, label="\$ x \$ - гарм.", histtype="step", linewidth=2, color="red")
        hist(zh, bins=[-2.0:0.05:2.0;], density=true, alpha=0.5, label="\$ z \$ - гарм.", histtype="step", linewidth=2, color="blue")
    end;
    xlabel("мкм")
    ylabel("Плотность вероятности")
    title("Распределение по координатам");
    legend()

    subplot(223)
    hist2D(vx, vz, bins=50, range=[[-0.3, 0.3], [-0.3, 0.3]],cmap="plasma", rasterized=true);
    xlabel("\$ v_x \$, мкм/мкс")
    ylabel("\$ v_z \$, мкм/мкс")
    title("Распределение по скоростям");

    subplot(224)
    hist(vx, bins=[-0.3:0.01:0.3;], density=true, alpha=0.5, label="\$ v_x \$")
    hist(vz, bins=[-0.3:0.01:0.3;], density=true, alpha=0.5, label="\$ v_z \$")
    if harmonic
        hist(vxh, bins=[-0.3:0.005:0.3;], density=true, alpha=0.5, label="\$ v_x \$ - гарм.", histtype="step", linewidth=2, color="red")
        hist(vzh, bins=[-0.3:0.005:0.3;], density=true, alpha=0.5, label="\$ v_z \$ - гарм.", histtype="step", linewidth=2, color="blue")
    end;
    xlabel("мкм/мкс")
    ylabel("Плотность вероятности")
    title("Распределение по скоростям");
    legend()

    tight_layout()
end;



function samples_visualise_eng(samples; samples_harmonic=[], harmonic=false, fontsize=16)
    x, y, z, vx, vy, vz = invert(samples);
    if harmonic
        xh, yh, zh, vxh, vyh, vzh = invert(samples_harmonic);
    end;

    figure(figsize=(10,10))
    subplot(221)
    hist2D(x, z, bins=50, range=[[-1.0, 1.0], [-2.5, 2.5]],cmap="plasma", rasterized=true);
    xlabel("x, μm", fontsize=fontsize)
    ylabel("z, μm", fontsize=fontsize)
    title("Coordinate distribution", fontsize=fontsize);

    subplot(222)
    hist(x, bins=[-1.5:0.02:1.5;], density=true, alpha=0.5, label="x")
    hist(z, bins=[-2.0:0.1:2.0;], density=true, alpha=0.5, label="z")
    if harmonic
        hist(xh, bins=[-1.5:0.01:1.5;], density=true, alpha=0.5, label="\$ x \$ - harmonic", histtype="step", linewidth=2, color="red")
        hist(zh, bins=[-2.0:0.05:2.0;], density=true, alpha=0.5, label="\$ z \$ - harmonic", histtype="step", linewidth=2, color="blue")
    end;
    xlabel("μm", fontsize=fontsize)
    ylabel("pdf", fontsize=fontsize)
    title("Velocity distribution", fontsize=fontsize);
    legend()

    subplot(223)
    hist2D(vx, vz, bins=50, range=[[-0.3, 0.3], [-0.3, 0.3]],cmap="plasma", rasterized=true);
    xlabel("\$ v_x \$, μm/μs", fontsize=fontsize)
    ylabel("\$ v_z \$, μm/μs", fontsize=fontsize)
    title("Velocity distribution", fontsize=fontsize);

    subplot(224)
    hist(vx, bins=[-0.3:0.01:0.3;], density=true, alpha=0.5, label="\$ v_x \$")
    hist(vz, bins=[-0.3:0.01:0.3;], density=true, alpha=0.5, label="\$ v_z \$")
    if harmonic
        hist(vxh, bins=[-0.3:0.005:0.3;], density=true, alpha=0.5, label="\$ v_x \$ - harmonic", histtype="step", linewidth=2, color="red")
        hist(vzh, bins=[-0.3:0.005:0.3;], density=true, alpha=0.5, label="\$ v_z \$ - harmonic", histtype="step", linewidth=2, color="blue")
    end;
    xlabel("μm/μs", fontsize=fontsize)
    ylabel("pdf", fontsize=fontsize)
    title("Velocity distribution", fontsize=fontsize);
    legend(fontsize=fontsize-2)

    tight_layout()
end;