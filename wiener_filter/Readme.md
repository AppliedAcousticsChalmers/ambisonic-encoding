# Multichannel Wiener Filter

This is our first attempt on a multichannel Wiener filter for suppressing noise in the microphone array signals. Using noise suppression is particular critical for the omni-based cardioid EMA.

The script estimates the spatial covariance of the background noise from the recording `background_noise_10s.wav` and then creates a Wiener filter based on that. The Wiener filter operates on the raw microphone array signals `voice_recording_noisy.wav` to create the file `voice_recording_denoised.wav`. These either of these signals can be rendered with the script `../render_ema_hybrid_to_ambisonics.m`.

It may have advantages to make it operate on the virtual cardioid signals that are created in the script `../render_ema_hybrid_to_ambisonics.m` instead of the raw microphone signals. We'll test that in the near future.