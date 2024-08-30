# esegui il file C
import os

import pandas as pd
import plotly.graph_objects as go

try:
    os.system(
        "rm /home/valiokei/Documenti/GitHub/nelder-mead/Crazyflie/Lib_YI/*CF.csv")
    os.system(
        "rm /home/valiokei/Documenti/GitHub/nelder-mead/Crazyflie/Lib_YI/main_CF_YI")
except:
    pass


os.system("/usr/bin/gcc -fdiagnostics-color=always -g /home/valiokei/Documenti/GitHub/nelder-mead/Crazyflie/Lib_YI/*CF_YI.c /home/valiokei/Documenti/GitHub/nelder-mead/Crazyflie/Lib_YI/*CF_YI.h -o /home/valiokei/Documenti/GitHub/nelder-mead/Crazyflie/Lib_YI/main_CF_YI -lm -I /home/valiokei/Documenti/GitHub/nelder-mead")

os.system(
    "cd /home/valiokei/Documenti/GitHub/nelder-mead/Crazyflie/Lib_YI/ && ./main_CF_YI")


# # Load CSV files into Pandas DataFrame
# df_Anchors = pd.read_csv('Anchors_CF_YI.csv', header=0)
# df_EstimatedTagPositions = pd.read_csv(
#     'EstimatedTagPositions_CF_YI.csv', header=0)
# df_TrueTagPositions = pd.read_csv('TrueTagPositions_CF_YI.csv', header=0)

# # add also the inputPoint File with the same format
# df_inputPoint = pd.read_csv('InputPoint_CF_YI.csv', header=0)


# # compute the euclidean distance between the estimated and true tag positions
# df_EstimatedTagPositions['Error'] = ((df_EstimatedTagPositions['Estimated_T_x'] - df_TrueTagPositions['True_T_x'])**2 + (
#     df_EstimatedTagPositions['Estimated_T_y'] - df_TrueTagPositions['True_T_y'])**2)**0.5

# # compute the mean error
# mean_error = df_EstimatedTagPositions['Error'].mean()
# print('Mean Error: ', mean_error)

# # Plot the data using Plotly
# fig = go.Figure()
# fig.add_annotation(
#     x=0.5, y=0.9, text=f"Mean Error: {mean_error:.2f}", showarrow=False)

# fig.add_trace(go.Scatter(
#     x=df_Anchors['AnchorX'], y=df_Anchors['AnchorY'], mode='markers', name='Anchors', marker=dict(color='green')))
# fig.add_trace(go.Scatter(
#     x=df_EstimatedTagPositions['Estimated_T_x'], y=df_EstimatedTagPositions['Estimated_T_y'], mode='markers', name='Estimated Tag Positions', marker=dict(color='blue')))
# fig.add_trace(go.Scatter(
#     x=df_TrueTagPositions["True_T_x"], y=df_TrueTagPositions["True_T_y"], mode='markers', name='True Tag Positions', marker=dict(color='red', symbol='circle-open')))
# fig.add_trace(go.Scatter(
#     x=df_inputPoint["InputPoint_x"], y=df_inputPoint["InputPoint_y"], mode='markers', name='Input Point', marker=dict(color='black', symbol='x')))


# fig.update_layout(title='Estimated Tag Positions',
#                   xaxis_title='X',
#                   yaxis_title='Y',
#                   width=900,
#                   height=500)  # Set the width and height to make the figure square

# # Add the mean error to the figure

# fig.show()
