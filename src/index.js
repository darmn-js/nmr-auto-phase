/**
 * Automatic Phase correction algorithm
 * @param {SD} spectraData - SD instance
 * @param {number} [phi0 = 0] - value
 * @param {number} [phi1 = 0] - value
 * @return {SD} returns the modified spectraData
 */
export default function autoPhaseCorrection(spectraData) {

  let nbPoints = spectraData.getNbPoints();
  let reData = spectraData.getYData(0);
  let imData = spectraData.getYData(1);
  let xData = spectraData.getXData(0);
  let magData = reData;//getMagnitudSpectrum(reData, imData);

  let ds = holoborodko(magData);
  let peaksDs = robustBaseLineRegionsDetection(ds);
  let peaksSp = robustBaseLineRegionsDetection(magData);
  let finalPeaks = new Array(reData.length);
  let xy = [];
  for (let i = 0; i < reData.length; i++) {
    finalPeaks[i] = peaksSp[i] | peaksSp[i];
    xy.push(xData[i]);
    xy.push(finalPeaks[i] * 1);
  }
  API.createData('mask', xy);
}

/**
 * Return a single array containing the magnitud of the complex array
 * @param {Array} re 
 * @param {Array} im 
 * @returns {Array} magnitud
 */
function getMagnitudSpectrum(re, im) {
  let mag = new Array(re.length);
  for (let i = 0; i < re.length; i++) {
    mag[i] = Math.sqrt(re[i] * re[i] + im[i] * im[i]);
  }
  return mag;
}

/**
 * Calculate the first derivative of the smoothed spectrum
 * @param {Array} s 
 * @returns {Array}  
 */
function holoborodko(s) {
  let dk = new Array(s.length);
  for (let i = 5; i < s.length - 5; i++) {
    dk[i] = (42 * (s[i + 1] - s[i - 1]) + 48 * (s[i + 2] - s[i - 2]) + 27 * (s[i + 3] + s[i - 3])
      + 8 * (s[i + 4] - s[i - 4]) + s[i + 5] - s[i - 5]) / 512;
  }
  //Fill the borders
  for (let i = 0; i < 5; i ++) {
      dk[i] = dk[5];
      dk[s.length - i - 1] = dk[s.length -  6]
  }
  
  return dk;
}

/**
 * Returns a binary mask, marking the zones containing peaks. 
 * @param {Array} s 
 * @returns {Array}
 */
function robustBaseLineRegionsDetection(s) {
  let mask = new Array(s.length);
  for (let i = 0; i < s.length; i++) {
    mask[i] = false;
  }

  //Recursivelly check for points greater than 3 times the sdt
  let change = true;
  while (change) {
    let res = stats(s, mask);
    let noiseLevel = 3 * res.std;
    console.log(noiseLevel)
    let mean = res.mean;
    change = false;
    for (let i = 0; i < s.length; i++) {
      if (Math.abs(s[i] - mean) > noiseLevel && !mask[i]) {
        change = true;
        mask[i] = true;
      }
    }
  }

  // Clean up mask by merging peaks blocks, separated by just a few points(4??).
  let count = 0;
  let prev = 0;
  for (let i = 0; i < s.length; i++) {
    if (!mask[i]) {
      count++;
    } else {
      if (count < 4) {
        for (let j = 0; j < count; j++) {
          mask[prev + j] = true;
        }
      }
      while (mask[++i] && i < s.length);
      prev = i;
      count = 0;
    }
  }

  return mask;
}

/**
 * Returns the mean of the spectrum
 * @param {Array} s 
 * @returns {number}
 */
function mean(s) {
  let sum = 0;
  for (let i = 0; i < s.length; i++) {
    sum += s[i];
  }
  return sum / s.length;
}

/**
 * Return the mean and the std of a given array, considering in the calculation only the points
 * where mask[i] == false
 * @param {Array} s 
 * @param {Array} mask 
 * @returns {object} {mean, std}
 */
function stats(s, mask) {
  let m = 0;
  let count = 0;
  for (let i = 0; i < s.length; i++) {
    if (!mask[i]) {
        m += s[i];
        count++;
    }
  }
  
  m /= count;
  let sum = 0;
 
  for (let i = 0; i < s.length; i++) {
    if (!mask[i]) {
      sum += Math.pow(s[i] - m, 2);
      count++;
    }
  }
  return { mean: m, std: Math.sqrt(sum / count) };
}
