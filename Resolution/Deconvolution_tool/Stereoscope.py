!pip install stereoscope
import stereoscope as st

sc_data = st.read("sc_data.h5ad")
spatial_data = st.read("spatial_data.h5ad")

# train the model
model = st.Stereoscope(sc_data, spatial_data)
model.fit()

# predict cell composition
results = model.predict()
