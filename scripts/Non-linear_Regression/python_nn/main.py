import matplotlib.pyplot as plt
import numpy as np
from keras.models import Sequential
from keras import optimizers
from keras.layers import Dense, Activation

# data points
X = np.arange(0.0, 5.0, 0.1, dtype='float32').reshape(-1,1)
y = 5 * np.power(X, 2) + np.power(np.random.randn(50).reshape(-1,1),3)

# model
model = Sequential()
model.add(Dense(50, activation='linear', input_dim=1))
model.add(Dense(30, activation='tanh', init='uniform'))
model.add(Dense(activation='linear', output_dim=1))

# training
sgd = optimizers.SGD(lr=0.001)
model.compile(loss='mse', optimizer=sgd)
model.fit(X, y, nb_epoch=1000)

# predictions
predictions = model.predict(X)

# plot
plt.scatter(X, y, edgecolors='g')
plt.plot(X, predictions,'r')
plt.legend(['Predicted Y', 'Actual Y'])
plt.show()