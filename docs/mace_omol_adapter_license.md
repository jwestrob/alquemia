# MACE interaction adapter attribution

`scripts/mace_omol_edges.py` adapts the operation order of
`RealAgnosticResidualNonLinearInteractionBlock.forward` from MACE0.3.16
(`mace/modules/blocks.py`), retaining the native modules and parameters while
batching neighbor edges. The source package's installed license notice follows.
This notice concerns the source implementation; checkpoint weights retain
their separately recorded license.

MIT License

Copyright (c) 2022 ACEsuit/mace

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
