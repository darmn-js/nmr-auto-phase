/**
 * Implementation of the algorithm for automatic phase correction: A robust, general automatic phase 
 * correction algorithm for high-resolution NMR data. 10.1002/mrc.4586
 * @param {object} spectraData 
 */
function autoPhaseCorrection(spectraData) {

  let nbPoints = spectraData.getNbPoints();
  let reData = spectraData.getYData(0);
  let imData = spectraData.getYData(1);
  let xData = spectraData.getXData(0);
  let magData = reData;//getMagnitudSpectrum(reData, imData);
  
  // TODO: It could be better to use the magnitud spectrum instead of the real data
  // for determining the peak regions
  //let magData = reData;//getMagnitudSpectrum(reData, imData);

  let ds = holoborodko(magData);
  let peaksDs = robustBaseLineRegionsDetection(ds);
  let peaksSp = robustBaseLineRegionsDetection(magData);
  let finalPeaks = new Array(nbPoints);
  let xy = [];
  for (let i = 0; i < nbPoints; i++) {
    finalPeaks[i] = peaksSp[i] & peaksSp[i];
    xy.push(xData[i]);
    xy.push(finalPeaks[i] * 1);
  }
  //API.createData('mask', xy);

  // Once the regions are detected, we auto phase each of them separately. 
  // TODO: This part can be put inside a function
  let i = -1;
  let x0 = 0;
  let res = [];
  while(i < nbPoints) {
      //phase first region
      let re = [];
      let im = [];
    
      //Look for the first 1 in the array
      while (!finalPeaks[++i] && i < nbPoints) {
        x0 = i;
      }
      
      //TODO: Add some extra points(0.1 ppm) at rigth and left sides of the region.
      while (finalPeaks[i] && i < nbPoints) {
        re.push(reData[i]);
        im.push(imData[i]);
        i++;
      }
      
      if (re.length > 0) {
        res.push(autoPhaseRegion(re, im, x0));
      }
  }
  
  // TODO: Still some corrections needed. In the paper they remove the outlayers interatively
  // until they can perform a regression witout bad points. Can someone help here?
  let reg = weightedLinearRegression(res.map(r => r.x0 / nbPoints), 
                                            res.map(r => r.ph0), 
                                            res.map(r => r.area));
                                            
  spectraData.phaseCorrection(reg[1] * Math.PI / 180, reg[0] * Math.PI / 180);

}

/**
 * Automatically determines the better phase correction of first order for the given spectrum.
 * Returns the first order phase correction to be applied in the region. It returns aswell the area of 
 * the spectrum and the index of the first point in the original spectrum.
 * @param {Array} re real part of the spectrum
 * @param {Array} im imaginary part of the spectrum
 * @param {number} x0 index of the first point in the original spectrum
 * @returns {object} {ph0, area, x0}
 */
function autoPhaseRegion(re, im, x0) {
  let start = -180;
  let stop = 180;
  let nSteps = 20;
  let maxSteps = 3;
  let bestAng = 0;
  while (maxSteps > 0) {
    let dAng = (stop - start) / (nSteps + 1);
    let minArea = Number.MAX_VALUE;
    bestAng = start;
    for (let i = start; i <= stop; i += dAng) {

      let phased = phaseCorrection(re, im, Math.PI * i / 180, 0);
      let negArea = 0;
      for (let j = 0; j < re.length; j++) {
        if (phased.re[j] < 0) {
          negArea += -phased.re[j];
        }
      }
      if (negArea < minArea) {
        minArea = negArea;
        bestAng = i;
      }
    }
    start = bestAng - dAng;
    stop = bestAng + dAng;
    maxSteps--;

  }
  
  // Calculate the area for the best angle
  let phased = phaseCorrection(re, im, Math.PI * bestAng / 180, 0);
  let area =0;
  for (let j = 0; j < re.length; j++) {
      area += phased.re[j];
  }
  //console.log(bestAng + " " + minArea);
  return { ph0: bestAng, area, x0 }
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
  for (let i = 0; i < 5; i++) {
    dk[i] = dk[5];
    dk[s.length - i - 1] = dk[s.length - 6]
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
  const SMALL = 64;
  for (let i = 0; i < s.length; i++) {
    if (!mask[i]) {
      count++;
    } else {
      if (count < SMALL) {
        for (let j = 0; j <= count; j++) {
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

/**
 * Perform a phase correction in the given data. Return new arrays containing the phased
 * spetrum
 * @param {Array} reData real data part
 * @param {Array} imData imaginary data part
 * @param {Number} phi0 in radians
 * @param {Number} phi1 in radians
 */
function phaseCorrection(reData, imData, phi0, phi1) {
  phi0 = Number.isFinite(phi0) ? phi0 : 0;
  phi1 = Number.isFinite(phi1) ? phi1 : 0;

  var nbPoints = reData.length;
  var re = new Array(nbPoints);
  var im = new Array(nbPoints);

  var delta = phi1 / nbPoints;
  var alpha = 2 * Math.pow(Math.sin(delta / 2), 2);
  var beta = Math.sin(delta);
  var cosTheta = Math.cos(phi0);
  var sinTheta = Math.sin(phi0);
  var cosThetaNew, sinThetaNew;

  var reTmp, imTmp;
  var index;
  for (var i = 0; i < nbPoints; i++) {
    index = nbPoints - i - 1;
    index = i;
    reTmp = reData[index] * cosTheta - imData[index] * sinTheta;
    imTmp = reData[index] * sinTheta + imData[index] * cosTheta;
    re[index] = reTmp;
    im[index] = imTmp;
    // calculate angles i+1 from i
    cosThetaNew = cosTheta - (alpha * cosTheta + beta * sinTheta);
    sinThetaNew = sinTheta - (alpha * sinTheta - beta * cosTheta);
    cosTheta = cosThetaNew;
    sinTheta = sinThetaNew;
  }

  return { re, im };
}

/**
 * Univariated weigthed linear regression. Only works for X having a single 
 * attribute. I've build the Newton's normal equations for this special case.
 * @param {Array} x independent data
 * @param {Array} y output
 * @param {Array} w weights
 * @returns {Array} [slope, intercept]
 */
function weightedLinearRegression(x, y, w) {
    let sxtw =0;
    let swx = 0;
    let sw = 0;
    let sxtwy = 0;
    let swy = 0;
    for (let i = 0; i< x.length; i++) {
        sxtw +=  x[i] * x[i] * w[i];
        swx += x[i] * w[i];
        sw += w[i];
        sxtwy += x[i] * w[i] * y[i];
        swy += w[i] * y[i];
    }
    
    /* Just to know what is the matrix system that we solve
     let Mx = [[sxtw, swx], [swx, sw]];
     let My = [[sxtwy], [swy]];
    */
    
    //Mx inverse 
    let detMx = sxtw * sw - swx * swx;
    let inMx = [[sw / detMx, -swx / detMx], [-swx / detMx, sxtw / detMx]];
    
    //let B = ml.MatrixLib.inverse(Mx).mmul(My);
    return [inMx[0][0] * sxtwy + inMx[0][1] * swy, 
             inMx[1][0] * sxtwy + inMx[1][1] * swy];
    
}
