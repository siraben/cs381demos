// PCA Iris Demo specific script

// UI Elements
const plotDivProjection = document.getElementById('plotProjection');
const plotDivMSE = document.getElementById('plotMSE');
const statusDiv = document.getElementById('status');

// Plotting Layouts
const PLOT_LAYOUT_PROJECTION = {
    title: 'Iris Data Projected onto PC1 & PC2',
    xaxis: { title: 'Principal Component 1' },
    yaxis: { title: 'Principal Component 2' },
    autosize: true,
    margin: { l: 40, r: 10, b: 40, t: 40 },
    legend: { title: { text: 'Species' } }
};
const PLOT_LAYOUT_MSE = {
    title: 'PCA Reconstruction MSE vs. Number of Components',
    xaxis: { title: 'Number of Principal Components (k)', tickvals: [1, 2, 3, 4] },
    yaxis: { title: 'Mean Squared Error (MSE)', autorange: true },
    autosize: true,
    margin: { l: 50, r: 10, b: 40, t: 40 }
};

// --- Plain JavaScript Math/Matrix Helpers ---

// Calculate column means of a matrix (array of arrays)
function calculateMeans(matrix) {
    if (!matrix || matrix.length === 0) return [];
    const nRows = matrix.length;
    const nCols = matrix[0].length;
    const means = Array(nCols).fill(0);
    for (let j = 0; j < nCols; j++) {
        for (let i = 0; i < nRows; i++) {
            means[j] += matrix[i][j];
        }
        means[j] /= nRows;
    }
    return means;
}

// Center data matrix by subtracting column means
function centerData(matrix, means) {
    if (!matrix || matrix.length === 0) return [];
    const nRows = matrix.length;
    const nCols = matrix[0].length;
    const centered = [];
    for (let i = 0; i < nRows; i++) {
        centered[i] = Array(nCols);
        for (let j = 0; j < nCols; j++) {
            centered[i][j] = matrix[i][j] - means[j];
        }
    }
    return centered;
}

// Calculate covariance matrix (sample covariance, divides by n-1)
function calculateCovariance(centeredMatrix) {
    if (!centeredMatrix || centeredMatrix.length < 2) return [];
    const nRows = centeredMatrix.length;
    const nCols = centeredMatrix[0].length;
    const cov = Array(nCols).fill(0).map(() => Array(nCols).fill(0));

    for (let j = 0; j < nCols; j++) {
        for (let k = j; k < nCols; k++) {
            let sum = 0;
            for (let i = 0; i < nRows; i++) {
                sum += centeredMatrix[i][j] * centeredMatrix[i][k];
            }
            cov[j][k] = sum / (nRows - 1);
            if (j !== k) {
                cov[k][j] = cov[j][k]; // Symmetric
            }
        }
    }
    return cov;
}

// Matrix multiplication C = A * B
function multiplyMatrices(A, B) {
    const rowsA = A.length, colsA = A[0].length;
    const rowsB = B.length, colsB = B[0].length;
    if (colsA !== rowsB) throw new Error("Matrix dimension mismatch for multiplication");

    const C = Array(rowsA).fill(0).map(() => Array(colsB).fill(0));
    for (let i = 0; i < rowsA; i++) {
        for (let j = 0; j < colsB; j++) {
            for (let k = 0; k < colsA; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Transpose matrix
function transpose(matrix) {
    if (!matrix || matrix.length === 0) return [];
    const rows = matrix.length, cols = matrix[0].length;
    const transposed = Array(cols).fill(0).map(() => Array(rows).fill(0));
    for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
            transposed[j][i] = matrix[i][j];
        }
    }
    return transposed;
}

// --- Jacobi Eigenvalue Algorithm for Symmetric Matrices ---
// Adapted from various sources (e.g., numeric.js, standard algorithms)
function jacobiEigenvalue(A, maxSweeps = 100, tolerance = 1e-10) {
    const n = A.length;
    let V = Array(n).fill(0).map((_, i) => Array(n).fill(0).map((__, j) => (i === j ? 1.0 : 0.0)));
    let D = A.map(row => [...row]); // Copy A, D will become diagonal
    let changed = true;

    for (let sweep = 0; sweep < maxSweeps && changed; sweep++) {
        changed = false;
        let sumOffDiagonalSq = 0;
        for (let p = 0; p < n; p++) {
            for (let q = p + 1; q < n; q++) {
                sumOffDiagonalSq += D[p][q] * D[p][q];
            }
        }

        if (Math.sqrt(sumOffDiagonalSq) < tolerance) break; // Converged

        for (let p = 0; p < n; p++) {
            for (let q = p + 1; q < n; q++) {
                if (Math.abs(D[p][q]) > tolerance / (n * n)) {
                    changed = true;
                    let t, c, s;
                    const app = D[p][p];
                    const aqq = D[q][q];
                    const apq = D[p][q];

                    if (Math.abs(apq) < tolerance * Math.abs(app - aqq)) {
                        t = apq / (aqq - app);
                    } else {
                        const theta = (aqq - app) / (2.0 * apq);
                        t = 1.0 / (Math.abs(theta) + Math.sqrt(theta * theta + 1.0));
                        if (theta < 0.0) t = -t;
                    }

                    c = 1.0 / Math.sqrt(t * t + 1.0);
                    s = t * c;

                    // Update D = R^T * D * R
                    D[p][p] = app + t * apq;
                    D[q][q] = aqq - t * apq;
                    D[p][q] = 0.0;
                    D[q][p] = 0.0;

                    for (let k = 0; k < n; k++) {
                        if (k !== p && k !== q) {
                            const akp = D[k][p];
                            const akq = D[k][q];
                            D[k][p] = c * akp + s * akq;
                            D[p][k] = D[k][p]; // Symmetric
                            D[k][q] = c * akq - s * akp;
                            D[q][k] = D[k][q]; // Symmetric
                        }
                    }

                    // Update V = V * R
                    for (let k = 0; k < n; k++) {
                        const vkp = V[k][p];
                        const vkq = V[k][q];
                        V[k][p] = c * vkp + s * vkq;
                        V[k][q] = c * vkq - s * vkp;
                    }
                }
            }
        }
    }

    const eigenvalues = Array(n);
    for (let i = 0; i < n; i++) {
        eigenvalues[i] = D[i][i];
    }

    // Sort eigenvalues (descending) and corresponding eigenvectors
    const eigenPairs = eigenvalues.map((value, index) => ({
        value: value,
        vector: V.map(row => row[index]) // Get column vector
    }));
    eigenPairs.sort((a, b) => b.value - a.value);

    return {
        eigenvalues: eigenPairs.map(p => p.value),
        eigenvectors: transpose(eigenPairs.map(p => p.vector)) // Store eigenvectors as columns
    };
}

// --- Main function using Plain JS ---
async function runPCAIrisDemo() {
    try {
        statusDiv.textContent = 'Loading data...';
        const response = await fetch('iris.json');
        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }
        const irisData = await response.json();
        statusDiv.textContent = 'Processing data...';

        const features = ['sepalLength', 'sepalWidth', 'petalLength', 'petalWidth'];
        const X_arr = irisData.map(d => features.map(f => d[f]));
        const y_arr = irisData.map(d => d.species);
        const speciesNames = [...new Set(y_arr)];
        const nSamples = X_arr.length;
        const nFeatures = features.length;

        // --- PCA Calculation using Plain JS ---
        statusDiv.textContent = 'Performing PCA...';

        // 1. Center the data
        const means = calculateMeans(X_arr);
        const X_centered = centerData(X_arr, means);

        // 2. Calculate Covariance Matrix
        const covMatrix = calculateCovariance(X_centered);

        // 3. Eigenvalue Decomposition using Jacobi Algorithm
        const { eigenvalues, eigenvectors } = jacobiEigenvalue(covMatrix);

        // 4. Project data onto first two principal components
        const W_2d = eigenvectors.map(row => [row[0], row[1]]); // Get first 2 columns
        const X_projected_2d_transposed = multiplyMatrices(transpose(W_2d), transpose(X_centered));
        const X_projected_arr = transpose(X_projected_2d_transposed);

        // --- Plot 1: 2D Projection (PC1 vs PC2) ---
        const projectionTraces = [];
        const colorMap = { 'setosa': 'red', 'versicolor': 'green', 'virginica': 'blue' };

        speciesNames.forEach(species => {
            const speciesIndices = y_arr.map((label, idx) => label === species ? idx : -1).filter(idx => idx !== -1);
            projectionTraces.push({
                x: speciesIndices.map(i => X_projected_arr[i][0]), // PC1
                y: speciesIndices.map(i => X_projected_arr[i][1]), // PC2
                mode: 'markers',
                type: 'scatter',
                name: species,
                marker: { color: colorMap[species], size: 7, opacity: 0.8 }
            });
        });

        Plotly.newPlot(plotDivProjection, projectionTraces, PLOT_LAYOUT_PROJECTION);
        statusDiv.textContent = 'Calculating Reconstruction MSE...';

        // --- Plot 2: Reconstruction MSE ---
        const mseValues = [];
        const componentsRange = [1, 2, 3, 4];

        for (let k of componentsRange) {
            // Select top k eigenvectors (columns)
            const W_k = eigenvectors.map(row => row.slice(0, k)); // Get first k columns
            const W_k_T = transpose(W_k);

            // Project onto k components: Z_k = W_k^T * X_centered^T
            const Z_k_T = multiplyMatrices(W_k_T, transpose(X_centered));

            // Reconstruct from k components: X_reconstructed^T = W_k * Z_k
            const X_reconstructed_centered_T = multiplyMatrices(W_k, Z_k_T);
            const X_reconstructed_centered = transpose(X_reconstructed_centered_T);

            // Calculate MSE
            let totalSquaredError = 0;
            for (let i = 0; i < nSamples; i++) {
                for (let j = 0; j < nFeatures; j++) {
                    const error = X_centered[i][j] - X_reconstructed_centered[i][j];
                    totalSquaredError += error * error;
                }
            }
            mseValues.push(totalSquaredError / (nSamples * nFeatures));
        }

        const mseTrace = {
            x: componentsRange,
            y: mseValues,
            mode: 'lines+markers',
            type: 'scatter',
            name: 'Reconstruction MSE'
        };

        // Adjust y-axis range for log scale if needed - REMOVED for linear scale
        /*
        if (mseValues[3] < 1e-15) {
           const minNonZeroMse = Math.max(1e-16, Math.min(...mseValues.filter(v => v > 1e-16)));
           PLOT_LAYOUT_MSE.yaxis.range = [Math.log10(minNonZeroMse) - 0.5, Math.log10(mseValues[0]) + 0.5];
        } else {
            PLOT_LAYOUT_MSE.yaxis.autorange = true;
        }
        */
       // Ensure autorange is set for linear scale (it's already there)
       PLOT_LAYOUT_MSE.yaxis.autorange = true;

        Plotly.newPlot(plotDivMSE, [mseTrace], PLOT_LAYOUT_MSE);

        statusDiv.textContent = 'Analysis complete.';

    } catch (error) {
        statusDiv.textContent = `Error: ${error.message}.`;
    }
}

// Run the demo
runPCAIrisDemo(); 