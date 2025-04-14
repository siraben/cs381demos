/**
 * plotUtils.js
 * Common utility functions for Plotly visualizations in the gradient descent demos.
 */

/**
 * Generates grid data for a 2D contour plot.
 * @param {function(number, number): number} f - The function f(w0, w1) to plot.
 * @param {Array<number>} range - Array [MIN, MAX] for the grid boundaries.
 * @param {number} step - The step size for the grid.
 * @returns {object} Object containing { w0vals, w1vals, zData }.
 */
function generateGridData(f, range, step) {
  const [MIN, MAX] = range;
  const w0vals = [];
  const w1vals = [];
  // Add a small epsilon to ensure MAX is included if step divides range evenly
  for (let w0 = MIN; w0 <= MAX + 1e-9; w0 += step) w0vals.push(w0);
  for (let w1 = MIN; w1 <= MAX + 1e-9; w1 += step) w1vals.push(w1);

  const zData = [];
  for (let i = 0; i < w1vals.length; i++) {
    const row = [];
    for (let j = 0; j < w0vals.length; j++) {
      row.push(f(w0vals[j], w1vals[i]));
    }
    zData.push(row);
  }
  return { w0vals, w1vals, zData };
}

/**
 * Creates a Plotly contour trace object.
 * @param {object} gridData - Object with { w0vals, w1vals, zData }.
 * @param {string} [colorscale='Viridis'] - Plotly colorscale name.
 * @param {number} [opacity=0.9] - Trace opacity.
 * @param {string} [name='Surface'] - Trace name.
 * @returns {object} Plotly trace configuration.
 */
function createContourTrace(gridData, colorscale = 'Viridis', opacity = 0.9, name = 'Surface') {
  return {
    x: gridData.w0vals,
    y: gridData.w1vals,
    z: gridData.zData,
    type: 'contour',
    colorscale: colorscale,
    contours: { coloring: 'heatmap', showlines: true },
    opacity: opacity,
    name: name
  };
}

/**
 * Creates a Plotly scatter trace object for the path.
 * @param {Array<number>} [x=[]] - Initial x coordinates.
 * @param {Array<number>} [y=[]] - Initial y coordinates.
 * @param {string} [color='red'] - Path color.
 * @param {number} [width=2] - Line width.
 * @param {number} [size=5] - Marker size.
 * @param {string} [name='Path'] - Trace name.
 * @returns {object} Plotly trace configuration.
 */
function createPathTrace(x = [], y = [], color = 'red', width = 2, size = 5, name = 'Path') {
  return {
    x: x,
    y: y,
    mode: 'lines+markers',
    line: { color: color, width: width },
    marker: { color: color, size: size },
    name: name
  };
}

/**
 * Creates a basic Plotly layout object for 2D plots.
 * @param {string} title - Plot title.
 * @param {Array<number>} range - Array [MIN, MAX] for axis ranges.
 * @param {string} [xaxisTitle='w0'] - X-axis title.
 * @param {string} [yaxisTitle='w1'] - Y-axis title.
 * @returns {object} Plotly layout configuration.
 */
function createLayout(title, range, xaxisTitle = 'w0', yaxisTitle = 'w1') {
  return {
    title: title,
    xaxis: { title: xaxisTitle, range: range, scaleratio: 1.0 }, // Maintain aspect ratio
    yaxis: { title: yaxisTitle, range: range, scaleratio: 1.0 }
  };
}

/**
 * Updates the text content of an element to display a slider's value.
 * @param {string} sliderId - The ID of the input range slider element.
 * @param {string} valueId - The ID of the element to display the value.
 * @param {boolean} [isInt=false] - If true, parse value as integer, else float with 2 decimals.
 * @param {function} [callback=null] - Optional callback function to execute after update.
 */
function setupSlider(sliderId, valueId, isInt = false, callback = null) {
  const slider = document.getElementById(sliderId);
  const valueDisplay = document.getElementById(valueId);
  if (!slider || !valueDisplay) {
      console.error(`Slider or value display not found for ${sliderId}/${valueId}`);
      return;
  }

  const updateValue = () => {
    const value = isInt ? parseInt(slider.value) : parseFloat(slider.value).toFixed(2);
    valueDisplay.textContent = value;
    if (callback) {
        callback(value);
    }
  };

  slider.addEventListener('input', updateValue);
  // Initial update
  updateValue();
}

/**
 * Initializes a Plotly plot.
 * @param {string} plotId - The ID of the div element for the plot.
 * @param {Array<object>} traces - An array of Plotly trace objects.
 * @param {object} layout - Plotly layout configuration.
 */
function plotInitial(plotId, traces, layout) {
  Plotly.newPlot(plotId, traces, layout);
}

/**
 * Updates the path data (x, y) of a specific trace in a Plotly plot.
 * @param {string} plotId - The ID of the plot div.
 * @param {Array<number>} newX - The new array of x coordinates.
 * @param {Array<number>} newY - The new array of y coordinates.
 * @param {number} [traceIndex=1] - The 0-based index of the trace to update (usually path is 1st or 2nd).
 */
function updatePath(plotId, newX, newY, traceIndex = 1) {
  Plotly.restyle(plotId, { x: [newX], y: [newY] }, [traceIndex]);
}

/**
 * Returns a random number within the plot range.
 * @param {Array<number>} range - [MIN, MAX]
 * @returns {number}
 */
function getRandomStart(range) {
    const [MIN, MAX] = range;
    return MIN + (MAX - MIN) * Math.random();
} 