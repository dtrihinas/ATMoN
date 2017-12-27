# ATMoN

ATMoN is an open-source framework developed to dynamically adjust the temporal granularity at which graph metrics are computed based on runtime knowledge captured by a low-cost probabilistic learning model capable of approximating both the metric stream evolution and the runtime volatility of the graph topology structure. This results in computationally offloading graph processing engines and eases the communication overhead in edge monitored networks over an unprecedented wealth of data.

Getting Started
---------------
The current ATMoN prototype is developed in R and requires the underlying open-source [igraph engine](http://igraph.org/r/).

How ATMoN works
--------------
coming soon...


Graph Metrics Available
-----------------------
- diameter
- effective diameter
- pagerank
- avg pagerank
- node outDegree
- avg outDegree
- betweeness centrality
- lobby index
- number of components
- giant component size
- max clique
- degree distribution
- number of nodes
- number of edges
- assortativity
- transitivity



Codebase Contributors
----------

- [Luis F. Chiroque](https://github.com/luisfo)
- [Demetris Trihinas](https://github.com/dtrihinas)

License and Disclaimer
----------
The framework is open-sourced under the Apache 2.0 License base. The codebase of the framework is maintained by the authors for academic research and is therefore provided "as is".

Reference
---------
When using the framework please use the following reference to cite our work:
TBD


Data
-------
The framework has been tested with several diverse and real-world datasets which can be acquired via direct downloads or after submitting a request. To obtain a list of datasets please send us an email (trihinas{at}cs.ucy.ac.cy)


