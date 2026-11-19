import Image from "next/image";
import Link from "next/link";
import { navGroups } from "./navigation";
import { assetPath } from "../site-path";

export default function DocsOverview() {
  return (
    <main className="docsMain docsOverview">
      <div className="docsArticle wide">
        <section className="overviewHero">
          <div className="overviewHeroCopy">
            <span className="docEyebrow">Hercules documentation</span>
            <h1>Hercules</h1>
            <p className="docLead">
              Hercules is a parallel, octree-based finite-element code for three-dimensional
              seismic wave propagation. Its adaptive mesh, distributed MPI execution, and
              GPU-accelerated solvers support simulations from regional wave propagation to
              local site and soil response.
            </p>
            <p className="overviewIntroduction">
              The physics extend beyond linear elasticity. Hercules can simulate nonlinear
              elastoplastic and hysteretic soils, including soil liquefaction with
              pressure-dependent constitutive models and nonlinear u-p coupling. It also
              supports saturated poroelastic media, acoustic domains, attenuation and
              absorbing layers, plus acoustic–elastic, poroelastic–elastic, nonlinear
              solid–pore-pressure, and acoustic–poroelastic coupled interfaces.
            </p>
          </div>
          <figure className="overviewHeroMark">
            <Image
              src={assetPath("/hercules-octree-mesh-mark.png?v=graded-centered-256")}
              alt="Hercules emblem formed from adaptive octree cells and a seismic waveform"
              width={512}
              height={512}
              priority
              unoptimized
            />
            <figcaption>Adaptive octree finite elements</figcaption>
          </figure>
        </section>
        <section className="docSection overviewCitation" aria-labelledby="cite-hercules-gpu">
          <h2 id="cite-hercules-gpu">Citing Hercules GPU</h2>
          <p>If you use Hercules GPU in your research, please cite this published workshop paper.</p>
          <p>
            Wenyang Zhang and Yinzhi Wang (2026).{" "}
            <cite>
              <a href="https://doi.ieeecomputersociety.org/10.1109/SCW72768.2026.00233">
                Modernizing Hercules for GPU-Accelerated Earthquake Ground Motion Simulation on NVIDIA GPUs
              </a>
            </cite>.{" "}
            <em>SC26 Workshops</em>, High Performance Computing for Environmental and Earth
            Sciences (HPC4EES). DOI 10.1109/SCW72768.2026.00233.
          </p>
        </section>
        <div className="overviewSummary">
          <p>This documentation covers the complete simulation workflow:</p>
          <ul>
            <li>building Hercules for an MPI and CUDA environment;</li>
            <li>preparing a case with one <code>parameters.in</code> file;</li>
            <li>running mesh generation and the solver on local or HPC systems;</li>
            <li>selecting elastic, nonlinear/plastic, liquefaction, acoustic, poroelastic, and coupled models;</li>
            <li>configuring outputs and performance paths; and</li>
            <li>validating numerical and physical results.</li>
          </ul>
        </div>
        <h2 className="overviewHeading">Contents</h2>
        <div className="docsOverviewGrid">
          {navGroups.map((group) => (
            <section key={group.label}>
              <h2>{group.label}</h2>
              {group.items.map((item) => (
                <Link href={`/docs/${item.slug}`} key={item.slug}>
                  <strong>{item.title}</strong>
                  <span>{item.description}</span>
                </Link>
              ))}
            </section>
          ))}
        </div>
      </div>
    </main>
  );
}
