import Link from "next/link";
import { navGroups } from "./navigation";

export default function DocsOverview() {
  return (
    <main className="docsMain docsOverview">
      <div className="docsArticle wide">
        <span className="docEyebrow">Hercules documentation</span>
        <h1>Hercules</h1>
        <p className="docLead">
          Hercules is a parallel, octree-based finite-element code for three-dimensional
          seismic wave propagation. It combines adaptive mesh generation, distributed
          execution, and GPU-accelerated solvers for linear, nonlinear elastoplastic,
          and coupled-physics simulations.
        </p>
        <div className="overviewSummary">
          <p>This documentation covers the complete simulation workflow:</p>
          <ul>
            <li>building Hercules for an MPI and CUDA environment;</li>
            <li>preparing a case with one <code>parameters.in</code> file;</li>
            <li>running mesh generation and the solver on local or HPC systems;</li>
            <li>selecting linear, nonlinear/plastic, and coupled physical models;</li>
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
