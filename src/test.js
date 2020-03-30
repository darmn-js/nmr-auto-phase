const autoPhaseCorrection = require('./index');

let spectrum = null;
if (get('zipData')) {
    let jcamp = get('zipData') //get("jcamp");
    let row = await processZip(jcamp);
    //console.log(row)
    let data = brukerConverter.convertFolder({
            acqus: DataObject.resurrect(row[0].acqus),
            "fid": row[0].content
        },
        {xy:true}
    );
    
    spectrum = new sd.NMR(data);
} else if(get("jcamp")) {
    spectrum = sd.NMR.fromJcamp(get("jcamp"));
}
console.log("GRPDLY " + spectrum.getParam("$GRPDLY"));

//spectrum.zeroFilling(spectrum.getNbPoints()*8)
//console.log(spectrum.getParam("$DECIM")+" " + spectrum.getParam("$DSPFVS"))
//spectrum.digitalFilter({nbPoints: -70})
spectrum.fourierTransform();
//Digital filtet in the fourier space
spectrum.phaseCorrection(0, spectrum.getParam("$GRPDLY") * (2 * Math.PI))//67.9842376708984*(2*Math.PI));
autoPhaseCorrection(spectrum);
API.createData("spectrum", spectrum.sd);

