import matplotlib.pyplot as plt
import numpy as np
from keras.models import Sequential
from keras import optimizers
from keras.layers import Dense
from math import sin

# Simulate Data
simulated_data = np.arange(-10.0, 10.0, 0.1, dtype='float32')
y = list(sin(x) / x for x in simulated_data)
simulated_data = simulated_data.reshape(-1, 1)
y = np.asarray(y).reshape(-1, 1)

# Define the model
model = Sequential()
model.add(Dense(10, activation='linear', input_dim=1))
model.add(Dense(10, activation='tanh', kernel_initializer='uniform'))
model.add(Dense(activation='linear', units=1))

# Train the Neural Network
sgd = optimizers.SGD(lr=0.01, nesterov=False)
model.compile(loss='mse', optimizer=sgd)
model.fit(simulated_data, y, epochs=30000)

# Make predictions
prediction = model.predict(simulated_data)

# Visualize
plt.figure(1)
plt.grid(True)
plt.scatter(simulated_data, y, edgecolors='g', linewidth=0.25)
plt.plot(simulated_data, prediction, 'r')
plt.legend(['Predicted', 'Actual'])
plt.xlabel('x')
plt.ylabel('y')
plt.title(r"NN with SGD Optimizer - Prediction of $\frac{sin(x)}{x}$ ")
# plt.savefig("figures/NN_SGD.pdf")
plt.show()
