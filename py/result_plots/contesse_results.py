# Copyright 2021 Adobe Research. All rights reserved.
# To view a copy of the license, visit LICENSE.md.

import streamlit as st
import pandas as pd
import altair as alt

st.title('ConTesse performance visualization')

ff = st.sidebar.checkbox("Front-faces only", value=True)

if ff:
	df = pd.read_csv('1017FF_all.csv')
	st.subheader('Front-faces only')
else:
	df = pd.read_csv('0927_all.csv')
	st.subheader('Front-faces and back-faces')


df = df.astype({'wso': 'string'})

df2 = pd.read_csv('0927faces_wso_stats.csv')
#df2.set_index('Input')

d1 = { df2['Input'][i]: int(df2['Control faces'][i]) for i in df2.index }
df['Quad faces']=list(df['name'].map(d1))

d2 = { df2['Input'][i]: int(df2['Input faces'][i]) for i in df2.index }
df['Lvl1 tris']=list(df['name'].map(d2))

#st.write(df)


font_size = st.sidebar.slider("Label font size",value=15,max_value=30)
title_font_size = st.sidebar.slider("Title font size",value=20,max_value=30)

c = alt.Chart(df).mark_circle().encode(
			x=alt.X('Lvl1 tris', title="Input triangles"),
			y=alt.Y('output_faces', title="Output triangles"),
			color=alt.Color('wso', title="Subdiv. levels"),
			tooltip=['name','cam','Lvl1 tris','output_faces'])

c=c.properties(title='Output size')
c=c.configure_axis(labelFontSize=font_size, titleFontSize=font_size)
c=c.configure_title(fontSize=title_font_size)
c=c.configure_legend(titleFontSize=font_size,labelFontSize=font_size) 



st.altair_chart(c,use_container_width=True)

c = alt.Chart(df).mark_circle().encode(
			x=alt.X('Lvl1 tris', title="Input triangles"),
			y=alt.Y('runtime', title="Runtime (seconds)"),
			color=alt.Color('wso', title="Subdiv. levels"),
			tooltip=['name','cam','Lvl1 tris','runtime'])

c=c.properties(title='Computation time')
c=c.configure_axis(labelFontSize=font_size, titleFontSize=font_size)
c=c.configure_title(fontSize=title_font_size)
c=c.configure_legend(titleFontSize=font_size,labelFontSize=font_size) 

st.altair_chart(c,use_container_width=True)