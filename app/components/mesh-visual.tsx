const nodes = [
  [10, 18], [26, 10], [44, 21], [62, 12], [81, 20],
  [18, 39], [36, 32], [54, 42], [73, 34], [91, 45],
  [8, 62], [29, 55], [47, 67], [67, 58], [84, 70],
  [17, 84], [39, 78], [59, 89], [78, 82], [94, 91],
];

export function MeshVisual() {
  return (
    <div className="meshVisual" aria-label="Abstract adaptive simulation mesh">
      <div className="meshGlow meshGlowOne" />
      <div className="meshGlow meshGlowTwo" />
      <div className="meshGrid" />
      <div className="waveArc waveArcOne" />
      <div className="waveArc waveArcTwo" />
      <div className="waveArc waveArcThree" />
      {nodes.map(([left, top], index) => (
        <span
          className={`meshNode ${index % 5 === 0 ? "meshNodeHot" : ""}`}
          style={{ left: `${left}%`, top: `${top}%` }}
          key={`${left}-${top}`}
        />
      ))}
      <div className="visualLegend">
        <span><i className="legendPulse" /> Wavefront</span>
        <span><i className="legendNode" /> Adaptive node</span>
      </div>
      <div className="visualReadout">
        <span>ELASTIC · ACOUSTIC · POROELASTIC</span>
        <strong>GPU-resident solve</strong>
      </div>
    </div>
  );
}
