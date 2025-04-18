// PCA 3D Reduction Demo specific script

// UI Elements
const plotDiv3d = document.getElementById('plot3d');
const plotDiv2d = document.getElementById('plot2d');
const numPointsSlider = document.getElementById('numPointsSlider');
const numPointsValue = document.getElementById('numPointsValue');
const spreadSlider = document.getElementById('spreadSlider');
const spreadValue = document.getElementById('spreadValue');
const resetBtn = document.getElementById('resetBtn');
const varianceInfo = document.getElementById('varianceInfo');

// Plotting Configuration
const PLOT_RANGE = [-8, 8]; // Range for all axes
const PLOT_LAYOUT_3D = {
    title: 'Original 3D Data & Principal Components',
    autosize: true,
    margin: { l: 0, r: 0, b: 0, t: 40 }, // Adjust margins
    scene: {
        xaxis: { title: 'x1', range: PLOT_RANGE },
        yaxis: { title: 'x2', range: PLOT_RANGE },
        zaxis: { title: 'x3', range: PLOT_RANGE },
        aspectratio: { x: 1, y: 1, z: 1 }
    }
};
const PLOT_LAYOUT_2D = {
    title: 'Data Projected onto PC1 & PC2',
    xaxis: { title: 'Principal Component 1', range: PLOT_RANGE, scaleratio: 1.0 },
    yaxis: { title: 'Principal Component 2', range: PLOT_RANGE, scaleratio: 1.0 },
    autosize: true,
    margin: { l: 40, r: 10, b: 40, t: 40 }
};

// --- Data Generation (3D) ---
function generate3DData(numPoints, spread) {
    const data = { x: [], y: [], z: [] };

    // Generate random rotation matrix (using Euler angles for simplicity)
    const alpha = Math.random() * 2 * Math.PI;
    const beta = Math.random() * Math.PI; // Only need 0 to PI
    const gamma = Math.random() * 2 * Math.PI;
    const ca = Math.cos(alpha), sa = Math.sin(alpha);
    const cb = Math.cos(beta), sb = Math.sin(beta);
    const cg = Math.cos(gamma), sg = Math.sin(gamma);

    const R = [
        [ca*cb*cg - sa*sg, -ca*cb*sg - sa*cg, ca*sb],
        [sa*cb*cg + ca*sg, -sa*cb*sg + ca*cg, sa*sb],
        [-sb*cg, sb*sg, cb]
    ];

    // Generate random scaling factors for axes
    const scale1 = (Math.random() * 0.6 + 0.4) * spread; // Stretch factor 1 (0.4 to 1.0 * spread)
    const scale2 = (Math.random() * 0.5 + 0.2) * spread; // Stretch factor 2 (0.2 to 0.7 * spread)
    const scale3 = (Math.random() * 0.4 + 0.1) * spread; // Stretch factor 3 (0.1 to 0.5 * spread)
    const S = [[scale1, 0, 0], [0, scale2, 0], [0, 0, scale3]];

    // Combine transformation T = R * S (matrix multiplication)
    const T = [[0,0,0],[0,0,0],[0,0,0]];
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            for (let k = 0; k < 3; k++) {
                T[i][j] += R[i][k] * S[k][j];
            }
        }
    }

    for (let i = 0; i < numPoints; i++) {
        // Generate standard normal variables (Box-Muller or just Math.random approximations)
        let z1 = Math.sqrt(-2.0 * Math.log(Math.random())) * Math.cos(2.0 * Math.PI * Math.random());
        let z2 = Math.sqrt(-2.0 * Math.log(Math.random())) * Math.sin(2.0 * Math.PI * Math.random());
        let z3 = Math.sqrt(-2.0 * Math.log(Math.random())) * Math.cos(2.0 * Math.PI * Math.random());

        // Apply transformation T
        let x = T[0][0] * z1 + T[0][1] * z2 + T[0][2] * z3;
        let y = T[1][0] * z1 + T[1][1] * z2 + T[1][2] * z3;
        let z_ = T[2][0] * z1 + T[2][1] * z2 + T[2][2] * z3; // Renamed to avoid conflict

        data.x.push(x);
        data.y.push(y);
        data.z.push(z_);
    }
    return data;
}

// --- PCA Calculation (3D) ---
function calculateMean3D(data) {
    if (!data || data.x.length === 0) return [0, 0, 0];
    const sumX = data.x.reduce((a, b) => a + b, 0);
    const sumY = data.y.reduce((a, b) => a + b, 0);
    const sumZ = data.z.reduce((a, b) => a + b, 0);
    const n = data.x.length;
    return [sumX / n, sumY / n, sumZ / n];
}

function calculateCovariance3D(data, mean) {
    if (!data || data.x.length < 2) return [[0,0,0],[0,0,0],[0,0,0]];
    const n = data.x.length;
    const [mX, mY, mZ] = mean;
    let cXX = 0, cXY = 0, cXZ = 0, cYY = 0, cYZ = 0, cZZ = 0;

    for (let i = 0; i < n; i++) {
        const dX = data.x[i] - mX;
        const dY = data.y[i] - mY;
        const dZ = data.z[i] - mZ;
        cXX += dX * dX;
        cXY += dX * dY;
        cXZ += dX * dZ;
        cYY += dY * dY;
        cYZ += dY * dZ;
        cZZ += dZ * dZ;
    }

    const factor = 1 / (n - 1);
    return [
        [cXX * factor, cXY * factor, cXZ * factor],
        [cXY * factor, cYY * factor, cYZ * factor], // Symmetric
        [cXZ * factor, cYZ * factor, cZZ * factor]  // Symmetric
    ];
}

// --- Eigenvalue/Eigenvector Calculation for 3x3 Symmetric Matrix ---
// Note: This involves solving a cubic equation and can be numerically sensitive.
// Using a library is generally preferred, but implementing directly for demo.

// Helper: Solve cubic equation x^3 + ax^2 + bx + c = 0
function solveCubic(a, b, c) {
    const Q = (a * a - 3 * b) / 9;
    const R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;
    const R2 = R * R;
    const Q3 = Q * Q * Q;

    let roots = [];
    if (R2 < Q3) { // Three real roots
        const theta = Math.acos(R / Math.sqrt(Q3));
        const sqrtQ = Math.sqrt(Q);
        roots.push(-2 * sqrtQ * Math.cos(theta / 3) - a / 3);
        roots.push(-2 * sqrtQ * Math.cos((theta + 2 * Math.PI) / 3) - a / 3);
        roots.push(-2 * sqrtQ * Math.cos((theta - 2 * Math.PI) / 3) - a / 3);
    } else { // One real root (or multiple equal roots)
        const A = -Math.sign(R) * Math.pow(Math.abs(R) + Math.sqrt(R2 - Q3), 1/3);
        const B = (A === 0) ? 0 : Q / A;
        roots.push((A + B) - a / 3); // The real root
        // For simplicity in this demo, we'll assume distinct eigenvalues or handle degeneracy loosely.
        // If A+B is the only real root, the other two are complex conjugates.
        // If R2 == Q3, there are multiple roots, but we only find one here easily.
        // This might cause issues if eigenvalues are very close.
    }
    // Sort roots numerically (descending for eigenvalues)
    roots.sort((x, y) => y - x);
    return roots;
}

// Helper: Calculate eigenvector for a given eigenvalue lambda
// Solves (A - lambda*I)v = 0 using cross products for 3x3 symmetric
function getEigenvector(A, lambda) {
    const B = [
        [A[0][0] - lambda, A[0][1], A[0][2]],
        [A[1][0], A[1][1] - lambda, A[1][2]],
        [A[2][0], A[2][1], A[2][2] - lambda]
    ];

    // Calculate cross products of rows of B
    const r1 = B[0], r2 = B[1], r3 = B[2];
    const c1 = crossProduct(r2, r3);
    const c2 = crossProduct(r3, r1);
    const c3 = crossProduct(r1, r2);

    // Find the cross product with the largest magnitude
    const m1 = dotProduct(c1, c1);
    const m2 = dotProduct(c2, c2);
    const m3 = dotProduct(c3, c3);

    let eigenvector;
    if (m1 >= m2 && m1 >= m3) {
        eigenvector = c1;
    } else if (m2 >= m1 && m2 >= m3) {
        eigenvector = c2;
    } else {
        eigenvector = c3;
    }

    return normalizeVector3D(eigenvector);
}

// Vector helpers
function crossProduct(a, b) {
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]];
}
function dotProduct(a, b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
function normalizeVector3D(v) {
    const len = Math.sqrt(dotProduct(v, v));
    if (len < 1e-9) return [0, 0, 0];
    return [v[0] / len, v[1] / len, v[2] / len];
}

function calculateEigen3x3(A) {
    // Characteristic polynomial: -lambda^3 + tr(A)*lambda^2 - k*lambda + det(A) = 0
    // Or lambda^3 - tr(A)*lambda^2 + k*lambda - det(A) = 0
    const m11 = A[0][0], m12 = A[0][1], m13 = A[0][2];
    const m22 = A[1][1], m23 = A[1][2];
    const m33 = A[2][2];

    const trace = m11 + m22 + m33;
    const det = m11*(m22*m33 - m23*m23) - m12*(m12*m33 - m23*m13) + m13*(m12*m23 - m22*m13);
    // k = sum of principal minors
    const k = (m11*m22 - m12*m12) + (m11*m33 - m13*m13) + (m22*m33 - m23*m23);

    // Solve lambda^3 - trace*lambda^2 + k*lambda - det = 0
    const eigenvalues = solveCubic(-trace, k, -det);

    // Calculate eigenvectors for each eigenvalue
    const eigenvectors = eigenvalues.map(lambda => getEigenvector(A, lambda));

    // Create pairs and sort by eigenvalue descending
    const eigenPairs = eigenvalues.map((val, index) => ({ value: val, vector: eigenvectors[index] }));
    eigenPairs.sort((a, b) => b.value - a.value);

    return {
        eigenvalues: eigenPairs.map(p => p.value),
        eigenvectors: eigenPairs.map(p => p.vector)
    };
}

// --- Data Projection ---
function projectData(data, mean, eigenvectors) {
    const [mX, mY, mZ] = mean;
    const v1 = eigenvectors[0]; // PC1
    const v2 = eigenvectors[1]; // PC2
    const n = data.x.length;
    const projected = { x: [], y: [] }; // x maps to PC1, y maps to PC2

    for (let i = 0; i < n; i++) {
        // Center the data point
        const centered = [data.x[i] - mX, data.y[i] - mY, data.z[i] - mZ];
        // Project onto PC1 and PC2
        projected.x.push(dotProduct(centered, v1));
        projected.y.push(dotProduct(centered, v2));
    }
    return projected;
}

// --- Plotting ---
function createPCA3DTraces(data, mean, pcaResult) {
    const dataTrace = {
        x: data.x, y: data.y, z: data.z,
        mode: 'markers',
        type: 'scatter3d',
        marker: { size: 3, color: 'blue', opacity: 0.6 },
        name: 'Data Points'
    };

    const [mX, mY, mZ] = mean;
    const { eigenvalues, eigenvectors } = pcaResult;
    const scaleFactor = 3; // Make vectors more visible
    const vectorTraces = [];

    for (let i = 0; i < 3; i++) {
        const ev = eigenvectors[i];
        const lambda = eigenvalues[i];
        const len = Math.sqrt(Math.max(0, lambda)) * scaleFactor;
        vectorTraces.push({
            x: [mX, mX + ev[0] * len], y: [mY, mY + ev[1] * len], z: [mZ, mZ + ev[2] * len],
            mode: 'lines',
            type: 'scatter3d',
            line: { color: i === 0 ? 'red' : (i === 1 ? 'green' : 'black'), width: 5 },
            name: `PC ${i+1}`
        });
    }

    return [dataTrace, ...vectorTraces];
}

function createProjected2DTrace(projectedData) {
    return {
        x: projectedData.x,
        y: projectedData.y,
        mode: 'markers',
        type: 'scatter',
        marker: { size: 5, color: 'green', opacity: 0.7 },
        name: 'Projected Data (PC1 vs PC2)'
    };
}

// --- Main Update Function ---
function updatePlots() {
    const numPoints = parseInt(numPointsSlider.value);
    const spread = parseFloat(spreadSlider.value);

    const data = generate3DData(numPoints, spread);
    const mean = calculateMean3D(data);
    const covMatrix = calculateCovariance3D(data, mean);
    const pcaResult = calculateEigen3x3(covMatrix);
    const projectedData = projectData(data, mean, pcaResult.eigenvectors);

    const traces3D = createPCA3DTraces(data, mean, pcaResult);
    const traces2D = [createProjected2DTrace(projectedData)];

    Plotly.react(plotDiv3d, traces3D, PLOT_LAYOUT_3D);
    Plotly.react(plotDiv2d, traces2D, PLOT_LAYOUT_2D);

    // Update variance info
    const totalVariance = pcaResult.eigenvalues.reduce((a, b) => a + b, 0);
    if (totalVariance > 1e-9) {
        const var1 = (pcaResult.eigenvalues[0] / totalVariance * 100).toFixed(1);
        const var2 = (pcaResult.eigenvalues[1] / totalVariance * 100).toFixed(1);
        const varTotal = ((pcaResult.eigenvalues[0] + pcaResult.eigenvalues[1]) / totalVariance * 100).toFixed(1);
        varianceInfo.textContent = `Variance captured by PC1: ${var1}%, PC2: ${var2}%, Total (PC1+PC2): ${varTotal}%`;
    } else {
        varianceInfo.textContent = 'Variance captured by PC1: -, PC2: -, Total: -';
    }
}

// --- Initialization and Event Listeners ---
function initialize() {
    setupSlider('numPointsSlider', 'numPointsValue', true, updatePlots);
    setupSlider('spreadSlider', 'spreadValue', false, updatePlots);
    resetBtn.addEventListener('click', updatePlots);
    updatePlots(); // Initial plot
}

// Run initialization
initialize(); 