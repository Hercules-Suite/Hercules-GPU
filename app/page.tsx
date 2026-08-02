import Link from "next/link";
import { MeshVisual } from "./components/mesh-visual";
import { SiteHeader } from "./components/site-header";

const entryPoints = [
  {
    index: "01",
    label: "First simulation",
    title: "Start with a working case",
    description:
      "Configure the toolchain, build Hercules, and complete a small end-to-end run.",
    href: "/docs/getting-started",
  },
  {
    index: "02",
    label: "Production",
    title: "Scale across GPUs",
    description:
      "Map ranks, use CUDA-aware MPI, keep fields resident, and select scalable I/O.",
    href: "/docs/gpu-execution",
  },
  {
    index: "03",
    label: "Reference",
    title: "Configure every path",
    description:
      "Understand inputs, physics, runtime controls, observations, and diagnostics.",
    href: "/docs/configuration",
  },
];

const capabilities = [
  ["Adaptive octrees", "Resolve complex three-dimensional structure without a uniformly fine grid."],
  ["Multiple physics", "Advance elastic, acoustic, and poroelastic domains with coupled interfaces."],
  ["GPU-resident fields", "Minimize host transfers through sparse synchronization and device interpolation."],
  ["Scalable execution", "Coordinate mesh construction, CUDA-aware MPI, and distributed output at scale."],
];

export default function Home() {
  return (
    <>
      <SiteHeader />
      <main>
        <section className="hero">
          <div className="heroBackdrop" />
          <div className="heroInner">
            <div className="heroCopy">
              <div className="eyebrowPill">
                <span /> GPU-accelerated wave propagation
              </div>
              <h1>
                Earthquake simulation,
                <br />
                <em>from mesh to insight.</em>
              </h1>
              <p>
                Hercules is an adaptive, distributed finite-element framework for
                high-fidelity seismic wave simulation—built for modern GPU clusters
                and complex coupled physics.
              </p>
              <div className="heroActions">
                <Link className="button buttonPrimary" href="/docs/getting-started">
                  Get started <span aria-hidden="true">→</span>
                </Link>
                <Link className="button buttonSecondary" href="/docs/architecture">
                  Explore the architecture
                </Link>
              </div>
              <div className="heroMeta" aria-label="Hercules capabilities">
                <span><strong>3D</strong> adaptive FEM</span>
                <span><strong>GPU</strong> resident solver</span>
                <span><strong>MPI</strong> distributed scale</span>
              </div>
            </div>
            <div className="heroVisualWrap">
              <MeshVisual />
            </div>
          </div>
        </section>

        <section className="entrySection">
          <div className="sectionHeading">
            <div>
              <span className="sectionKicker">Choose your path</span>
              <h2>Documentation built around the work</h2>
            </div>
            <p>
              Go from a first local case to a reproducible production simulation
              without losing the physical and computational assumptions along the way.
            </p>
          </div>
          <div className="entryGrid">
            {entryPoints.map((item) => (
              <Link className="entryCard" href={item.href} key={item.index}>
                <div className="cardTopline">
                  <span>{item.index}</span>
                  <small>{item.label}</small>
                </div>
                <h3>{item.title}</h3>
                <p>{item.description}</p>
                <span className="cardArrow" aria-hidden="true">↗</span>
              </Link>
            ))}
          </div>
        </section>

        <section className="capabilitySection">
          <div className="capabilityIntro">
            <span className="sectionKicker light">One computational framework</span>
            <h2>Complex physics.<br />Clear execution paths.</h2>
            <p>
              The documentation connects physical models to the mesh, communication,
              synchronization, and output decisions that affect every result.
            </p>
            <Link href="/docs/physics">Explore physics models →</Link>
          </div>
          <div className="capabilityList">
            {capabilities.map(([title, description], index) => (
              <article key={title}>
                <span>0{index + 1}</span>
                <div>
                  <h3>{title}</h3>
                  <p>{description}</p>
                </div>
              </article>
            ))}
          </div>
        </section>

        <section className="quickstartBand">
          <div>
            <span className="sectionKicker">A reproducible beginning</span>
            <h2>Start with the launcher you will use at scale.</h2>
            <p>
              The unified workflow builds or imports the mesh, selects GPU ranks, and
              transfers mesh leaves in memory within one allocation.
            </p>
          </div>
          <div className="terminalCard">
            <div className="terminalBar">
              <span><i /> <i /> <i /></span>
              <small>run.sh</small>
              <button aria-label="Code sample shown below">shell</button>
            </div>
            <pre><code><span># Configure the allocation</span>{"\n"}export BASE_IN=/scratch/case{"\n"}export PHASE_A_NP=64{"\n"}export PHASE_B_NP=8{"\n\n"}<b>./run_phaseA_phaseB.sh</b></code></pre>
          </div>
        </section>
      </main>
      <footer className="siteFooter">
        <div>
          <span className="brandMark small" aria-hidden="true">H</span>
          <p>Hercules Documentation<br /><small>Adaptive finite-element wave simulation</small></p>
        </div>
        <nav aria-label="Footer navigation">
          <Link href="/docs/getting-started">Documentation</Link>
          <Link href="/docs/validation">Benchmarks</Link>
          <a href="https://github.com/zwyll/GPU_new">Source ↗</a>
        </nav>
      </footer>
    </>
  );
}
