// PCA Demo specific script

// UI Elements
const plotDiv = document.getElementById('plot');
const numPointsSlider = document.getElementById('numPointsSlider');
const numPointsValue = document.getElementById('numPointsValue');
const spreadSlider = document.getElementById('spreadSlider');
const spreadValue = document.getElementById('spreadValue');
const resetBtn = document.getElementById('resetBtn');

// Plotting Configuration
const PLOT_X_RANGE = [-6, 6];
const PLOT_Y_RANGE = [-6, 6];
const PLOT_LAYOUT = createLayout('PCA Visualization', PLOT_X_RANGE, 'x1', 'x2'); // Assuming createLayout is in plotUtils.js

let currentData = { x: [], y: [] };

// --- Data Generation ---
/**
 * Generates 2D normally distributed data points with a random covariance structure.
 * @param {number} numPoints - Number of points to generate.
 * @param {number} spread - Controls the overall standard deviation.
 */
function generateData(numPoints, spread) {
    const data = { x: [], y: [] };
    // Create a random covariance matrix (symmetric positive semi-definite)
    // To make it interesting, generate data skewed along a random direction
    const angle = Math.random() * Math.PI;
    const cosA = Math.cos(angle);
    const sinA = Math.sin(angle);
    const scale1 = (Math.random() * 0.5 + 0.5) * spread; // Stretch factor 1 (0.5 to 1.0 * spread)
    const scale2 = (Math.random() * 0.5 + 0.2) * spread; // Stretch factor 2 (0.2 to 0.7 * spread)

    // Rotation matrix R
    const R = [[cosA, -sinA], [sinA, cosA]];
    // Scaling matrix S
    const S = [[scale1, 0], [0, scale2]];
    // Transformation matrix T = R * S
    const T = [
        [R[0][0] * S[0][0], R[0][1] * S[1][1]],
        [R[1][0] * S[0][0], R[1][1] * S[1][1]]
    ];

    for (let i = 0; i < numPoints; i++) {
        // Generate standard normal variables (mean 0, std dev 1)
        // Box-Muller transform
        let u1 = Math.random();
        let u2 = Math.random();
        let z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
        let z2 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2);

        // Apply transformation T to z = [z1, z2]
        let x = T[0][0] * z1 + T[0][1] * z2;
        let y = T[1][0] * z1 + T[1][1] * z2;

        data.x.push(x);
        data.y.push(y);
    }
    currentData = data;
    return data;
}

// --- PCA Calculation ---
/**
 * Calculates the mean of the data.
 * @param {object} data - { x: Array<number>, y: Array<number> }
 * @returns {Array<number>} [meanX, meanY]
 */
function calculateMean(data) {
    if (!data || data.x.length === 0) return [0, 0];
    const sumX = data.x.reduce((a, b) => a + b, 0);
    const sumY = data.y.reduce((a, b) => a + b, 0);
    return [sumX / data.x.length, sumY / data.y.length];
}

/**
 * Calculates the 2x2 covariance matrix of the data.
 * @param {object} data - { x: Array<number>, y: Array<number> }
 * @param {Array<number>} mean - [meanX, meanY]
 * @returns {Array<Array<number>>} 2x2 covariance matrix [[covXX, covXY], [covYX, covYY]]
 */
function calculateCovariance(data, mean) {
    if (!data || data.x.length < 2) return [[0, 0], [0, 0]]; // Need at least 2 points
    const n = data.x.length;
    const [meanX, meanY] = mean;
    let covXX = 0, covXY = 0, covYY = 0;

    for (let i = 0; i < n; i++) {
        const diffX = data.x[i] - meanX;
        const diffY = data.y[i] - meanY;
        covXX += diffX * diffX;
        covXY += diffX * diffY;
        covYY += diffY * diffY;
    }

    // Use n-1 for sample covariance
    covXX /= (n - 1);
    covXY /= (n - 1);
    covYY /= (n - 1);

    return [[covXX, covXY], [covXY, covYY]]; // covYX = covXY
}

/**
 * Calculates the eigenvectors and eigenvalues of a 2x2 symmetric matrix.
 * @param {Array<Array<number>>} matrix - [[a, b], [b, d]]
 * @returns {object} { eigenvalues: [l1, l2], eigenvectors: [[v1x, v1y], [v2x, v2y]] }
 * Eigenvalues/vectors are sorted by eigenvalue magnitude (descending).
 */
function calculateEigen2x2(matrix) {
    const a = matrix[0][0];
    const b = matrix[0][1];
    const d = matrix[1][1];

    // Quadratic equation for eigenvalues: lambda^2 - trace*lambda + det = 0
    const trace = a + d;
    const det = a * d - b * b;

    // Eigenvalues
    const discriminant = Math.sqrt(trace * trace - 4 * det);
    const lambda1 = (trace + discriminant) / 2;
    const lambda2 = (trace - discriminant) / 2;

    let ev1, ev2;

    // Eigenvectors: (A - lambda*I)v = 0
    // For lambda1:
    // (a - lambda1)x + by = 0
    // bx + (d - lambda1)y = 0
    // Can use either equation. Handle case where b is near zero.
    if (Math.abs(b) > 1e-6 || Math.abs(a - lambda1) > 1e-6) {
        ev1 = normalizeVector([-b, a - lambda1]);
    } else { // If first row is essentially [0, 0]
        ev1 = normalizeVector([d - lambda1, -b]);
    }

    // For lambda2:
    if (Math.abs(b) > 1e-6 || Math.abs(a - lambda2) > 1e-6) {
        ev2 = normalizeVector([-b, a - lambda2]);
    } else {
        ev2 = normalizeVector([d - lambda2, -b]);
    }

    // Ensure orthogonality numerically if needed (should be guaranteed for symmetric)

    // Sort by eigenvalue magnitude (descending)
    if (Math.abs(lambda1) >= Math.abs(lambda2)) {
        return { eigenvalues: [lambda1, lambda2], eigenvectors: [ev1, ev2] };
    } else {
        return { eigenvalues: [lambda2, lambda1], eigenvectors: [ev2, ev1] };
    }
}

/**
 * Normalizes a 2D vector (makes its length 1).
 * @param {Array<number>} v - [x, y]
 * @returns {Array<number>} Normalized vector [nx, ny]
 */
function normalizeVector(v) {
    const len = Math.sqrt(v[0] * v[0] + v[1] * v[1]);
    if (len < 1e-9) return [0, 0]; // Avoid division by zero
    return [v[0] / len, v[1] / len];
}

// --- Plotting ---
/**
 * Creates Plotly traces for PCA visualization.
 * @param {object} data - { x: Array<number>, y: Array<number> }
 * @param {Array<number>} mean - [meanX, meanY]
 * @param {object} pcaResult - { eigenvalues: [l1, l2], eigenvectors: [[v1x, v1y], [v2x, v2y]] }
 * @returns {Array<object>} Array of Plotly trace objects.
 */
function createPCATraces(data, mean, pcaResult) {
    const dataTrace = {
        x: data.x,
        y: data.y,
        mode: 'markers',
        type: 'scatter',
        marker: { size: 6, color: 'blue', opacity: 0.7 },
        name: 'Data Points'
    };

    const [meanX, meanY] = mean;
    const { eigenvalues, eigenvectors } = pcaResult;
    const [ev1, ev2] = eigenvectors;
    const [lambda1, lambda2] = eigenvalues;

    // Scale vectors by sqrt(eigenvalue) for visualization (represents std dev)
    const scaleFactor = 2; // Make vectors more visible
    const v1_scaled_len = Math.sqrt(Math.max(0, lambda1)) * scaleFactor;
    const v2_scaled_len = Math.sqrt(Math.max(0, lambda2)) * scaleFactor;

    // Create vectors starting from the mean
    const vectorTrace = {
        x: [meanX, meanX + ev1[0] * v1_scaled_len,
            meanX, meanX + ev2[0] * v2_scaled_len],
        y: [meanY, meanY + ev1[1] * v1_scaled_len,
            meanY, meanY + ev2[1] * v2_scaled_len],
        mode: 'lines',
        type: 'scatter',
        line: { color: 'red', width: 3 },
        name: 'Principal Components'
    };

    // Separate trace for mean point
    const meanTrace = {
        x: [meanX],
        y: [meanY],
        mode: 'markers',
        marker: { size: 8, color: 'red', symbol: 'cross' },
        name: 'Data Mean'
    };

    return [dataTrace, vectorTrace, meanTrace];
}

// --- Main Update Function ---
function updatePlot() {
    const numPoints = parseInt(numPointsSlider.value);
    const spread = parseFloat(spreadSlider.value);

    const data = generateData(numPoints, spread);
    const mean = calculateMean(data);
    const covMatrix = calculateCovariance(data, mean);
    const pcaResult = calculateEigen2x2(covMatrix);

    const traces = createPCATraces(data, mean, pcaResult);

    Plotly.react(plotDiv, traces, PLOT_LAYOUT); // Use react for potentially faster updates
}

// --- Initialization and Event Listeners ---
function initialize() {
    // Setup sliders using utility function (assuming plotUtils.js is loaded)
    setupSlider('numPointsSlider', 'numPointsValue', true, updatePlot);
    setupSlider('spreadSlider', 'spreadValue', false, updatePlot);

    // Button listener
    resetBtn.addEventListener('click', updatePlot);

    // Initial plot
    updatePlot();
}

// Run initialization when the script loads
initialize(); 