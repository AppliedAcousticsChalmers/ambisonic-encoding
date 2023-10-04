# Ambisonic Encoding of Signals From Higher-Order Microphone Arrays

NB: 
* You will need 'git LFS' if you want to checkout all of the large data files. Install it with `git lfs install`. If you are still experiencing problems, use `git lfs pull` instead of `git pull`.
* See the branch 'octave' for a version that works in Octave.

## Spherical Microphone Arrays With a Spherical Baffle

The MATLAB script `render_sma_to_ambisonics.m` demonstrates how to compute ambisonic signals from the signals that are captured by the microphones of a classical spherical microphone array with a rigid spherical baffle. 

The underlying concept was presented in the literature in many locations in different variants. The present implementation uses definitions of the involved quantities such that the resulting ambisonic signals ares compatible with software tools like [SPARTA](https://leomccormack.github.io/sparta-site/) and the [IEM Plugin Suite](https://plugins.iem.at/). The mathematical formulation is summarized in

> J. Ahrens, "Ambisonic Encoding of Signals From Spherical Microphone Arrays," Technical note v. 1, Chalmers University of Technology, Aug. 2022 [[pdf]](https://arxiv.org/pdf/2211.00583.pdf).

If you execute the script then the microphone signals from [this recording](https://youtu.be/qcqeygqjxZ4?t=31) will be encoded into ambisonics. Additionally, a binaural preview will be computed, which is a binaural rendering of those ambisonic signals while assuming that the listener is looking straight. Note that the result of this computation is already included in the repository (`out_ambisonics.wav` and `out_sma_binaural.wav`) (Actually, by default, `out_ambisonics.wav` contains the encoding from the EMA presented here below. That one is much higher order, and I like it even more). You'll need to download the employed HRIRs from [here](https://zenodo.org/record/3928297/files/HRIR_L2702.sofa?download=1) and store them in the subfolder `resources` (The MATLAB script is going to do that automatically for you.) as well as the SOFA MATLAB API from [here](https://sourceforge.net/projects/sofacoustics/) for being able to compute the binaural preview yourself. The ambisonic encoding works either way.

The repository comprises also the [Reaper](https://www.reaper.fm/) project `binaural_rendering.rpp` that reads the ambisonic signals and renders them binaurally in realtime. This allows for using head tracking if you happen to have a tracker available. The Reaper project requires the [IEM Plugin Suite](https://plugins.iem.at/) to be installed. Note that an evaluation version of Reaper is sufficient for this. If you use a [Supperware headtracker](https://supperware.co.uk/headtracker-overview), you’ll need the *Bridgehead*, which you can obtain from the manufacturer’s website. Go to `Bridge settings` and set the `Profile` to `IEM Scene Rotator`. 

In Reaper, 

- Go to `Preferences -> Control/OSC/web -> Add`
- Set the `Control surface mode` to `OSC (Open Sound Control)` 
- Set `Mode` to `Local port [receive only]`
- Set `Local listen port`: `8000`

Note that the example implementation of the rendering that we provide does not comprise any form of equalization, which is usually useful to mitigate the effects of spherical harmonic order truncation and spatial aliasing. Thomas Deppisch provides a MATLAB implementation of eMagLS on [his GitHub repository](https://github.com/thomasdeppisch/eMagLS) that does exactly this. 

## Equatorial Microphone Arrays With a Spherical Baffle

The MATLAB script `render_ema_to_ambisonics.m` demonstrates how to compute ambisonic signals from the signals that are captured by the microphones of an equatorial microphone array (EMA). The underlying concept was originally presented in

> J. Ahrens, H. Helmholz, D. L. Alon, S. V. Amengual Garí, “Spherical Harmonic Decomposition of a Sound Field Based on Observations Along the Equator of a Rigid Spherical Scatterer” in J. Acoust. Soc. Am. 150(2), 2021 [[pdf]](http://www.soundfieldsynthesis.org/wp-content/uploads/pubs/Ahrens_etal_JASA2021.pdf).

The MATLAB script uses a reformulation of the equatorial array solution that is again compatible with software tools like [SPARTA](https://leomccormack.github.io/sparta-site/) and the [IEM Plugin Suite](https://plugins.iem.at/). The details are described in 

> J. Ahrens, "Ambisonic Encoding of Signals From Equatorial Microphone Arrays," Technical note v. 1, Chalmers University of Technology, Aug. 2022 [[pdf]](https://arxiv.org/pdf/2211.00584.pdf).

If you execute the scripts then the microphone signals from [this recording](https://youtu.be/95qDd13pVVY?t=58) will be encoded into 7th-order ambisonics and stored in the file `out_ambisonics.wav`, which is the same file that is mentioned above (it is getting overwritten). Start Reaper after computing the signals so that it loads the updated ones. The binaural preview will be stored in `out_ema_binaural.wav`.

eMagLS for EMAs is also available in [Tommi's repository](https://github.com/thomasdeppisch/eMagLS). 

We thank Reality Labs for funding the initial work on the EMA concept.

## License
The content of this repository is licensed under the terms of the MIT license. Please consult the file [LICENSE](LICENSE) for more information about this license.