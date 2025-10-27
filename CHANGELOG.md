# CHANGELOG


## v2.1.5 (2025-10-27)

### Bug Fixes

* fix: wrong nodata set with downloaded band data, caused issues in
clipping to AOI ([`77981b9`](https://github.com/ollinevalainen/satellitetools/commit/77981b93c2bb52292ea58bbf485c8ed7041a1cf4))

### Unknown

* Merge branch 'master' into develop ([`f724a10`](https://github.com/ollinevalainen/satellitetools/commit/f724a109d384c2237125e7e72bc4ccd6561ba18b))


## v2.1.4 (2025-10-23)

### Bug Fixes

* fix: ensure request params have default values if None given ([`1bd15cd`](https://github.com/ollinevalainen/satellitetools/commit/1bd15cd24c39e972359757d388766339b4582d83))


## v2.1.3 (2025-09-26)

### Bug Fixes

* fix: updated packages due to dependabot security notifications ([`04a7699`](https://github.com/ollinevalainen/satellitetools/commit/04a769933bbd2fe2703344a48af2b1072ac909a4))

* fix: ensure that defaults set if None given to request parameters ([`f54367d`](https://github.com/ollinevalainen/satellitetools/commit/f54367d3fffa2a9fea9fa55edc59514e70df30be))


## v2.1.2 (2025-02-13)

### Bug Fixes

* fix: multipolygon for aws and multipolygon tests ([`6fa9365`](https://github.com/ollinevalainen/satellitetools/commit/6fa9365159a376c11c65e695e739cfc8491552ff))

* fix: changed aws resampling approach ([`2e63964`](https://github.com/ollinevalainen/satellitetools/commit/2e6396472104d1550e82ab332e54a6eeae5f6836))

* fix: added missing 10 m band (B8) ([`1650558`](https://github.com/ollinevalainen/satellitetools/commit/1650558eb44a64eedbfa0dcffe362193cd3a441d))

* fix: crs property key in different pystac verisons ([`016d17f`](https://github.com/ollinevalainen/satellitetools/commit/016d17fc280fcd71841516244cb226a5085af134))

* fix: handling of missing band data in GEE images ([`fae9469`](https://github.com/ollinevalainen/satellitetools/commit/fae9469ce5e3ea2fc43d8bcee3d813c008fa96d5))

* fix: failing aws test to not use multiprocessing due to some random
occuring rasterio error ([`9d83eb7`](https://github.com/ollinevalainen/satellitetools/commit/9d83eb74cc6973afabb2e46e941cfdaf45a3e269))

### Documentation

* docs: added citing information ([`e563771`](https://github.com/ollinevalainen/satellitetools/commit/e56377145356ab5180ee4479f67254180f0c705e))

* docs: Updated todo and versions history to README ([`abf3415`](https://github.com/ollinevalainen/satellitetools/commit/abf34150cdea483f96d67d9b80051d6afbd7f681))

### Unknown

* Merge pull request #5 from juliusvira/allow-multipolygon-master

fix: Allow multipolygon master ([`15f17f0`](https://github.com/ollinevalainen/satellitetools/commit/15f17f0c357767bd029d7fd32420a70af63b8509))

* Merge branch 'ollinevalainen:master' into allow-multipolygon-master ([`92cecd3`](https://github.com/ollinevalainen/satellitetools/commit/92cecd38355bf563fa9dae25ecd4bb334e9a85e2))

* change list of geometries to just geometry ([`cae638f`](https://github.com/ollinevalainen/satellitetools/commit/cae638f9d5f29b2265ef953a2d3a135956115e45))

* multipolygon in the new gee ([`4a6f76b`](https://github.com/ollinevalainen/satellitetools/commit/4a6f76be61c20376afbaa159f17b6d27605497e4))


## v2.1.1 (2025-02-11)

### Bug Fixes

* fix: changed aws resampling approach ([`2204825`](https://github.com/ollinevalainen/satellitetools/commit/220482510b65616d78da53f5ef2dba2697e24210))

* fix: added missing 10 m band (B8) ([`9e731f5`](https://github.com/ollinevalainen/satellitetools/commit/9e731f5d8714ecdf83dd8ba321b92ffa364d776c))

* fix: crs property key in different pystac verisons ([`b7acdba`](https://github.com/ollinevalainen/satellitetools/commit/b7acdbadcf503c377275e1bc967c2eecbf12a5df))

* fix: handling of missing band data in GEE images ([`7af22b1`](https://github.com/ollinevalainen/satellitetools/commit/7af22b1b9430cbdd7a9e6b041487d758348695b5))

### Documentation

* docs: added citing information ([`d18e813`](https://github.com/ollinevalainen/satellitetools/commit/d18e813945da16fbde1607d36d437e8b74941a75))

* docs: Updated todo and versions history to README ([`fcf8125`](https://github.com/ollinevalainen/satellitetools/commit/fcf8125ef63b46d222631dee3a59371696590303))


## v2.1.0 (2025-01-08)

### Bug Fixes

* fix: main to master typo in workflow ([`e4ca6ab`](https://github.com/ollinevalainen/satellitetools/commit/e4ca6abcc2cc9386587c33ff31bbb436a1056a58))

* fix: missin init in test_wrappers, authentication success ([`9295cd9`](https://github.com/ollinevalainen/satellitetools/commit/9295cd9977df8f74d54c1f48527a9afb6d1b74a5))

* fix: issues with decoding base64 ([`0add63c`](https://github.com/ollinevalainen/satellitetools/commit/0add63cfe4ff2965c0e56a229fd9901b07af5246))

* fix: trying file again ([`9679e69`](https://github.com/ollinevalainen/satellitetools/commit/9679e69373baf8ce1389c09aa632c9585b3b5518))

* fix: bug in secret name ([`3b226fe`](https://github.com/ollinevalainen/satellitetools/commit/3b226fe4ba4fbe0ca20ce06fa6377bef4035408a))

* fix: debugging ([`189f73e`](https://github.com/ollinevalainen/satellitetools/commit/189f73e5e8b4b19c4fec87cac4e2549a7ff4aac8))

* fix: issues with the json file ([`3ca1c37`](https://github.com/ollinevalainen/satellitetools/commit/3ca1c3794eeca213c10bb4f7cefdeedae2464591))

* fix: trying key json from base64 secret ([`28db542`](https://github.com/ollinevalainen/satellitetools/commit/28db542cd55a0cd7aca64c17da910d1e6dc7df82))

* fix: added handling of GEE key file or data ([`17749ac`](https://github.com/ollinevalainen/satellitetools/commit/17749ac6494bde80c66da73cff9091fdfc460289))

* fix: bug in coverage ([`4579cd2`](https://github.com/ollinevalainen/satellitetools/commit/4579cd268c8de322f49839e72eb25cf48bfe182d))

* fix: debugging workflow ([`7864220`](https://github.com/ollinevalainen/satellitetools/commit/786422078bb6ea4c390ab368e12718ed562701b3))

* fix: debugging workflow ([`0168f2e`](https://github.com/ollinevalainen/satellitetools/commit/0168f2e3fac96c7b2feeaad3fd68baceb74107ab))

* fix: changed env location in workflow ([`43b53e1`](https://github.com/ollinevalainen/satellitetools/commit/43b53e1a80ee8a63e40aa757f5a3e2da520e7d32))

* fix: GEE service account credentials for workflow ([`ee020ff`](https://github.com/ollinevalainen/satellitetools/commit/ee020ff186d5266f40af61ac937d0aa2f74d03a8))

* fix: added pytest-cov ([`2132bed`](https://github.com/ollinevalainen/satellitetools/commit/2132bed2137e41d1ea1bc7c3f091db261e88f35b))

* fix: abandoned nox approach ([`5c42ebd`](https://github.com/ollinevalainen/satellitetools/commit/5c42ebd84f26fd97450be7f0304fcfd84c0e4b17))

* fix: bug in workflow definition ([`e6cf93b`](https://github.com/ollinevalainen/satellitetools/commit/e6cf93b13dea9313a2b1ccaac796082edc02614a))

### Build System

* build: added CI workflow ([`6e57f81`](https://github.com/ollinevalainen/satellitetools/commit/6e57f814b9af440708dd115111c274026fff4ee7))

### Features

* feat: set-up CD workflow ([`ce25914`](https://github.com/ollinevalainen/satellitetools/commit/ce25914bd7e0cf7ec8c13d0b12281e2114339c46))

### Unknown

* fix: ([`d2f732d`](https://github.com/ollinevalainen/satellitetools/commit/d2f732d9009249d069dd55683130691fc2e2d93c))

* Changed one print to logger.info ([`27fc1e9`](https://github.com/ollinevalainen/satellitetools/commit/27fc1e9b4a7dcb0256a456a866b0f5f0ba9d2002))


## v2.0.0 (2024-12-30)

### Unknown

* Prepared for publishing and releasing ([`25e9707`](https://github.com/ollinevalainen/satellitetools/commit/25e9707e3e02627d6fffd82d80e028cc8011c4a2))

* Requirements for RTD ([`c84b074`](https://github.com/ollinevalainen/satellitetools/commit/c84b074a02fb704f268e9b55b86af41ff6afb8ed))

* Added RTD config ([`961ffa1`](https://github.com/ollinevalainen/satellitetools/commit/961ffa1c377cfa7d99eed77edbd69e4045932d72))

* Updated README ([`0692980`](https://github.com/ollinevalainen/satellitetools/commit/0692980b49935d683f9bc31406d557b255107bfb))

* Updated README ([`2c1edfa`](https://github.com/ollinevalainen/satellitetools/commit/2c1edfac5f0c1fefcf93837d63009e321c3b71a2))

* From print to logging ([`4e67510`](https://github.com/ollinevalainen/satellitetools/commit/4e6751095da2bad7a3d7860d2bf4390523b83186))

* Updated README and docs ([`c16a3cf`](https://github.com/ollinevalainen/satellitetools/commit/c16a3cf8edec9ae7e316d52489148e412318b8c3))

* Updated todo list ([`565a73e`](https://github.com/ollinevalainen/satellitetools/commit/565a73e0e7e29c1e873cece3a8509c78ba3f9e9a))

* Created documentation with sphinx ([`35385b9`](https://github.com/ollinevalainen/satellitetools/commit/35385b9e18526861c080acc7c0f193d4d83fa6f8))

* Updated to-do list ([`6996574`](https://github.com/ollinevalainen/satellitetools/commit/69965744ab624d4435240a24fd1888b442c3b5a0))

* Fixd bug when x or y dim len 1 ([`3f037b6`](https://github.com/ollinevalainen/satellitetools/commit/3f037b6c01b5aac1c069c5156487beef5bca26ac))

* Contains filter to avoid aoi not fully contained within tile ([`b6601ca`](https://github.com/ollinevalainen/satellitetools/commit/b6601cae6f0c2b12a25cf87ec1c55524f9713390))

* Updated README ([`626da1f`](https://github.com/ollinevalainen/satellitetools/commit/626da1fe63098e615e78e1b3fbcb08588382bf13))

* Run tests for python 3.12 with nox ([`c74e926`](https://github.com/ollinevalainen/satellitetools/commit/c74e92651265f63f1e05720eba0c769833b8f52f))

* Added nox and run tests for 3.10 and 3.11 ([`b4dae2f`](https://github.com/ollinevalainen/satellitetools/commit/b4dae2f38bea9faa5a41bda3bea6990d372a9e0e))

* Updated README ([`5e3caf1`](https://github.com/ollinevalainen/satellitetools/commit/5e3caf104fc29b4c5f08f49455955b3abd49338f))

* Merge branch 'master' into develop ([`4e1709d`](https://github.com/ollinevalainen/satellitetools/commit/4e1709d06da91007fa73298bd199805aa5b6a522))

* Fixed breaking change in str Enum in python > 3.10, tested in 3.10 and 3.11 ([`ee7a95f`](https://github.com/ollinevalainen/satellitetools/commit/ee7a95f0991060a945fdb875c4b4e2d6e93fd7de))

* Updated README ([`a1e25d2`](https://github.com/ollinevalainen/satellitetools/commit/a1e25d24ec69e0949e5494dc571c8e4342314ab6))

* Updated README ([`72788b4`](https://github.com/ollinevalainen/satellitetools/commit/72788b40e6cb0e0732d78519fd4157b4abd1be1c))

* Fixed wrapper pipelines and improved printing ([`07b2476`](https://github.com/ollinevalainen/satellitetools/commit/07b247696c711ca27807eeffcc70d7be14f3e436))

* acquisition_time to utc aware ([`90cd13f`](https://github.com/ollinevalainen/satellitetools/commit/90cd13f26570794b2c170c480b08d80d0461bb7e))

* Merge branch 'develop' ([`d15e34d`](https://github.com/ollinevalainen/satellitetools/commit/d15e34d9f16dec4ba69e97b97d948aefac73ce5d))

* Added validation and type conversion for bands ([`56599d1`](https://github.com/ollinevalainen/satellitetools/commit/56599d1cd13f669f2fd4ff21cd1f0a0cf5bc2399))

* Merge branch 'master' of github.com:ollinevalainen/satellitetools ([`3484fb3`](https://github.com/ollinevalainen/satellitetools/commit/3484fb306b906a32a4408b7e8067d6bf02ae05cb))

* Update README.md ([`75f5e84`](https://github.com/ollinevalainen/satellitetools/commit/75f5e845efb96cfc5234c09133ee6f753f00e3ac))

* Removed unnecessary files and packages ([`cc99176`](https://github.com/ollinevalainen/satellitetools/commit/cc9917654061aa3c962e9855a96e64497f8c064e))

* satelliettools to version 2.0.0 ([`4546ac0`](https://github.com/ollinevalainen/satellitetools/commit/4546ac0a67907d5aa717af9dcf35e9958b637b34))

* Updated README ([`af7b670`](https://github.com/ollinevalainen/satellitetools/commit/af7b67097622b0a31092fec29f4e7207e1e0433d))

* Merge branch 'develop' ([`c0d362b`](https://github.com/ollinevalainen/satellitetools/commit/c0d362b2d7736651e1a3ab27d6cc13043ddb69c2))

* Fixed a bug in aoi_nan_percentage reported by Johan Sparf ([`551e7ae`](https://github.com/ollinevalainen/satellitetools/commit/551e7aef49dd68e7cb658fa32001226eebe07cc1))

* Updated version constraints due to security issues reported by
dependabot ([`fe7f876`](https://github.com/ollinevalainen/satellitetools/commit/fe7f8763b496fac1fbe59ddcee138e14efdf37cb))

* Updated README ([`55aac3e`](https://github.com/ollinevalainen/satellitetools/commit/55aac3ef5939dd0475a3cf5d47d19cdd359ae442))

* Cleaning code, refactored s2 item filtering ([`0fd1b27`](https://github.com/ollinevalainen/satellitetools/commit/0fd1b27ddbebb529c286709d50f667e166b1f866))

* Search only to GEE, biophys bug, sentinel2 tests ([`d5ed39c`](https://github.com/ollinevalainen/satellitetools/commit/d5ed39c5c52fd79d6fd25533ead02f83500bd3fd))

* Relaxed dependencies ([`433fad3`](https://github.com/ollinevalainen/satellitetools/commit/433fad3a38d4e64aea35a6409dabb91fa2a170f9))

* Check for missing SCL data, multiprocessing for wrapper, commenting ([`4e13077`](https://github.com/ollinevalainen/satellitetools/commit/4e13077daf02538f260154401f7a8c6a7e824194))

* Split long queries to 6 month intervals ([`3c3256a`](https://github.com/ollinevalainen/satellitetools/commit/3c3256a3d0dfafca98addaffce0c346b1e72b8bc))

* Increased default request limit ([`0780ef8`](https://github.com/ollinevalainen/satellitetools/commit/0780ef8f12a0844ded325d5004efa632c9e7bd06))

* Added data consistency check ([`5860335`](https://github.com/ollinevalainen/satellitetools/commit/5860335c57b15838a9bb88866b436334de873859))

* Added test for 2020 ([`61e7f1d`](https://github.com/ollinevalainen/satellitetools/commit/61e7f1d25ba36de22ba6cc622cfabf970e4b93b0))

* Updated README ([`a1adb37`](https://github.com/ollinevalainen/satellitetools/commit/a1adb375173f34a8a78c8d44937eea99f61d054b))

* Added test for wrappers ([`3fe4d37`](https://github.com/ollinevalainen/satellitetools/commit/3fe4d370940d9f64b7726df9f73b4804e3c8d05f))

* Updated README, changed and added methods, AWS_COG to AWS ([`7399a82`](https://github.com/ollinevalainen/satellitetools/commit/7399a82cb7f2e4902cc5f893600db42c25819684))

* Tested and fixed wrappers ([`1436c26`](https://github.com/ollinevalainen/satellitetools/commit/1436c26dd220b577dd535f262bd2a920501e15d7))

* Handling 2022 data collection for AWS ([`18fa042`](https://github.com/ollinevalainen/satellitetools/commit/18fa0429db9ffccf58bac40f91c9267fc101e60d))

* Major restructuring of code base. Improved AWS approach. Added tests. ([`410ad79`](https://github.com/ollinevalainen/satellitetools/commit/410ad7933fee9143af8c4a1a6fb85cffad054a64))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`82a3554`](https://github.com/ollinevalainen/satellitetools/commit/82a35541be1e4f64016d953fc84656a4130c22f2))

* Update README.md ([`803b54a`](https://github.com/ollinevalainen/satellitetools/commit/803b54a647dde6cf92f873bd852a9f039e044c60))

* Removed ee.Initialize() due to changed EarthEngine authentication. User
should run ee.Initialize themselves ([`ffcec1a`](https://github.com/ollinevalainen/satellitetools/commit/ffcec1acf171c1a4196bbb276e09aa192efec206))

* Removed old uncertainty function ([`23609f6`](https://github.com/ollinevalainen/satellitetools/commit/23609f650ef1a4bbb2881ec7e42b989e98bd33eb))

* Merge branch 'develop-poetry' into develop ([`9839992`](https://github.com/ollinevalainen/satellitetools/commit/9839992eb4037f1f4f35f2d51c5cffd7ac07a9a3))

* Changed to new Earth Search collection sentinel-2-c1-l2a ([`892520c`](https://github.com/ollinevalainen/satellitetools/commit/892520cb84e5fdc49aa829851cda138943eb4028))

* Function to get Copernicus DEM value for lat lon ([`f58fefa`](https://github.com/ollinevalainen/satellitetools/commit/f58fefaecb11d0dd14da5eb9ef91e403cde54829))

* Changed from concat to df from list of dicts ([`d3b35ad`](https://github.com/ollinevalainen/satellitetools/commit/d3b35adcaab80acb15a1ad9fbb8fe03a3415c744))

* Removed deprecated append for pandas frames ([`c75e8b4`](https://github.com/ollinevalainen/satellitetools/commit/c75e8b4c77d5dc794a6b9daa53d2684fe00e3558))

* Upgraded geopandas ([`592a46d`](https://github.com/ollinevalainen/satellitetools/commit/592a46d70b18d6bea37af888e0bc5a564892ccbd))

* Upgraded pandas ([`b1c2b9e`](https://github.com/ollinevalainen/satellitetools/commit/b1c2b9eac1d75d90c56777b43db89586934402c7))

* Update README.md

Installation command ([`730f223`](https://github.com/ollinevalainen/satellitetools/commit/730f22352639d4f0e9c0098c286cc1d9c35219d0))

* Changed from unmaintained sat-search to pystac-client ([`19453fd`](https://github.com/ollinevalainen/satellitetools/commit/19453fdc06bf732ae6927bb105c84bf3c77303e1))

* Merge branch 'develop' into develop-poetry ([`f4575ae`](https://github.com/ollinevalainen/satellitetools/commit/f4575ae329de9dbcc560ff72c6b596feb9245166))

* Fixed imports and updated examples ([`0468f6c`](https://github.com/ollinevalainen/satellitetools/commit/0468f6c465c65face69bae6a7a1aba7a8eb60465))

* Changed __init__ files ([`7fcccb0`](https://github.com/ollinevalainen/satellitetools/commit/7fcccb0568c26654eb8eb5f1867e017a61becf30))

* Passing pre-commits ([`787a478`](https://github.com/ollinevalainen/satellitetools/commit/787a4788b98b7604cfb3e50f9f806846afea2085))

* Conformed to flake8 partially. Not passing at the moment ([`702689c`](https://github.com/ollinevalainen/satellitetools/commit/702689c25633202caf458bc54d556f6a8d68e903))

* Updated requiremetns without hashes ([`98a9dca`](https://github.com/ollinevalainen/satellitetools/commit/98a9dca6a72937aa72823ce7821a2c99ddab0227))

* Created requirements.txt ([`3bdc223`](https://github.com/ollinevalainen/satellitetools/commit/3bdc2232f4a60dd801c33db2b1629988c1a26675))

* Installed pre-commit hooks ([`ee74816`](https://github.com/ollinevalainen/satellitetools/commit/ee74816cfe6fb0677d9c0aeee3f5b2cdb6bcde92))

* Added dev dependencies ([`239aff7`](https://github.com/ollinevalainen/satellitetools/commit/239aff7493a96afad9aaf46b911f2cf5b432c788))

* Non-conflicting dependencies found ([`6320e5f`](https://github.com/ollinevalainen/satellitetools/commit/6320e5f805a0e5661b557eb730cffd87ece4e31f))

* Trying to solve dependency conflicts, not solved ([`a3405ed`](https://github.com/ollinevalainen/satellitetools/commit/a3405edd912d5a4d98755164d9ca1ba3a77ebc49))

* Removed examples ([`aa1f19b`](https://github.com/ollinevalainen/satellitetools/commit/aa1f19b377533d11b239f332b86727afc5c8f3f3))

* Conformed to new folder structure ([`b9a1832`](https://github.com/ollinevalainen/satellitetools/commit/b9a1832b1914cade0152d8d366d35d4bb62ee485))

* Initialized poetry project ([`6b2f571`](https://github.com/ollinevalainen/satellitetools/commit/6b2f57112232e413fbac47f71e26f826cfeb3746))

* Updated uncertainty for small areas ([`b55503c`](https://github.com/ollinevalainen/satellitetools/commit/b55503cc5f9875f3aabc71d73e7127cf0a9aba04))

* Updated ucertainty formulation from SE to STD ([`c84402b`](https://github.com/ollinevalainen/satellitetools/commit/c84402bddfff3d67b10dd8316658258d5769adf2))

* Fxied bug in band name mappings ([`dab47c9`](https://github.com/ollinevalainen/satellitetools/commit/dab47c956b22dc53e271bbad976b50a306c9358a))

* Adapted AWS source to Earth serach v1, older v0 deprecated ([`926631c`](https://github.com/ollinevalainen/satellitetools/commit/926631cef179bfb996ae6c7f5bc11ab8ac6e7df5))

* Changed to S2_SR_HARMONIZED dataset due to new S2 process baseline ([`628ce87`](https://github.com/ollinevalainen/satellitetools/commit/628ce87e437f8d1eac4a33cd09713194fa50ae66))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`088c9b7`](https://github.com/ollinevalainen/satellitetools/commit/088c9b7e5bb37c326ac859045d7d627cd0c50e32))

* Added definition domain files for FCOVER ([`5a406f1`](https://github.com/ollinevalainen/satellitetools/commit/5a406f18cf18b2bb988f301504e44b7013449bef))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`54d2aad`](https://github.com/ollinevalainen/satellitetools/commit/54d2aadabdbba3bff9eb597f74b8831df4728f9a))

* Option to define scale for qi evaluation ([`09c9b2a`](https://github.com/ollinevalainen/satellitetools/commit/09c9b2a2b55d87c7292db164fbc25221735cffb3))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`8b4cdc2`](https://github.com/ollinevalainen/satellitetools/commit/8b4cdc2f2e85341929ca3e5ee89019803fe56b2b))

* Update README.md ([`8bb622a`](https://github.com/ollinevalainen/satellitetools/commit/8bb622a50744745c5a6a08c544e1f884e2f30e22))

* Fixed bug when all available data is inconsistent (bands are not equal
length) ([`a0bfcfe`](https://github.com/ollinevalainen/satellitetools/commit/a0bfcfe794cab5f82c5ef4ed6e86cc3af44d613f))


## v1.0.0 (2022-01-11)

### Unknown

* Changed deprecated to_wkt() to wkt attribute ([`a14e030`](https://github.com/ollinevalainen/satellitetools/commit/a14e030a82c1c9d1a86272d808b9f7d3259b1c07))

* Fixed bug when adjusted n gets too small ([`f0345ac`](https://github.com/ollinevalainen/satellitetools/commit/f0345ac3af29d789edf86ce87e1c7a5d2c8350ed))

* Updated GPP model usage ([`c36cb68`](https://github.com/ollinevalainen/satellitetools/commit/c36cb686bb6d87d5268890706d71f2144144b4c0))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`28756df`](https://github.com/ollinevalainen/satellitetools/commit/28756dfdb45b5bc7e021b68cddde92b8d8eea2f5))

* Update README.md ([`0da69d0`](https://github.com/ollinevalainen/satellitetools/commit/0da69d0c7d693776da64569373e56123bac3408d))

* Improved example ([`5eb40f0`](https://github.com/ollinevalainen/satellitetools/commit/5eb40f0fee169a399748aa7090c5d63e0925331f))

* Added new filter ([`88cce5b`](https://github.com/ollinevalainen/satellitetools/commit/88cce5b4d434fb2118dd7989f35694dbaa92ba27))

* Added code to remove inconsistent raster data ([`38dbaaf`](https://github.com/ollinevalainen/satellitetools/commit/38dbaaff149040726f209e85c70d0155d21f7473))

* Drop rows where SCL data is nan ([`bc72fd6`](https://github.com/ollinevalainen/satellitetools/commit/bc72fd6f1d2275f05acec51a5d240c13efe356e9))

* Added full projection information ([`1385a27`](https://github.com/ollinevalainen/satellitetools/commit/1385a270cf98e69e57a9ad75c0baf53822534c25))

* Merge branch 'develop-harmonize-gee-aws' into develop

Conflicts:
	aws_cog/aws_cog.py ([`9deb5f6`](https://github.com/ollinevalainen/satellitetools/commit/9deb5f607d6a8ade3a34a05eacf489bde67c7692))

* Added datasource to qi dataframes ([`f02e595`](https://github.com/ollinevalainen/satellitetools/commit/f02e5955520151cb371c0f460bf22c3e5428c541))

* Added utc=True to Date handling ([`280c96e`](https://github.com/ollinevalainen/satellitetools/commit/280c96ecbd5cd2510a564afca1be3b60f6663857))

* Added datasource as variable to gee ([`71cc389`](https://github.com/ollinevalainen/satellitetools/commit/71cc3897ffcada1a3d006d124f316222758a7af1))

* Added gsd param to gee+fixed aws pipeline issues ([`3163034`](https://github.com/ollinevalainen/satellitetools/commit/3163034bf19ffe16c32bb2033ae65d0a4eb96692))

* Removed bad import ([`80540da`](https://github.com/ollinevalainen/satellitetools/commit/80540da8f8d0eb2448f182d8c1930652c38a422e))

* Fixed masking and resampling issues ([`4b1c85e`](https://github.com/ollinevalainen/satellitetools/commit/4b1c85e41e6b1cd3d372dba360fd1fdfa4182945))

* Fixed bugs in update masking. Harmonizing gee-aws ([`4baae6d`](https://github.com/ollinevalainen/satellitetools/commit/4baae6dd92e0b0f8ce8b685c924ea6f579250d4b))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`9330e75`](https://github.com/ollinevalainen/satellitetools/commit/9330e75e183df1f633badcc696b1469ba80185b6))

* Fixed aws netcdf coords from pix corner to center ([`c9a68ab`](https://github.com/ollinevalainen/satellitetools/commit/c9a68abf1225721a5dabc006bad52448e2cb2135))

* Added clipping of nodata from raster borders ([`b8ba246`](https://github.com/ollinevalainen/satellitetools/commit/b8ba246ee57c6440ac5cd50100095552f1d17943))

* Increased window expansion from 50m to 80m ([`5fffd0c`](https://github.com/ollinevalainen/satellitetools/commit/5fffd0c24cedc6298cd01fa53933cc621891b98e))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`1281d0d`](https://github.com/ollinevalainen/satellitetools/commit/1281d0dc0a17a3db4222540bb2bc4e526b8fced2))

* Fixed n due to resampling. ([`a98e137`](https://github.com/ollinevalainen/satellitetools/commit/a98e1378aca6597657185f0c23da269b6ef8de23))

* Updated uncertainty computation ([`adafced`](https://github.com/ollinevalainen/satellitetools/commit/adafced7982dc504f648889ccf0bc11716b25815))

* Added confidence interval calculation ([`462a18e`](https://github.com/ollinevalainen/satellitetools/commit/462a18e98bce45b72005e892c189a2593fb4794a))

* Working on fix for uncertainties ([`2c909ab`](https://github.com/ollinevalainen/satellitetools/commit/2c909ab121bc50418dec47c31a05f8fae1dfcc98))

* Included previous uncert.code just in case ([`5c736a4`](https://github.com/ollinevalainen/satellitetools/commit/5c736a4b0eebc27ebffeaa94e911f78e67a1b16e))

* Added confidence level calculation ([`9cc046b`](https://github.com/ollinevalainen/satellitetools/commit/9cc046bed8d699597aeb9ead1a452e18c310534b))

* Fixed SCL gsd and aws_cog transformation errors ([`eb86abc`](https://github.com/ollinevalainen/satellitetools/commit/eb86abcc6abebebeacad248e999864f1377a9c71))

* did something ([`1983dc8`](https://github.com/ollinevalainen/satellitetools/commit/1983dc84f6a124a8cc2fa64527ecd33cd15e1047))

* Added uncertainty calculation for timeseries data ([`bdcc226`](https://github.com/ollinevalainen/satellitetools/commit/bdcc2265982e32229909767beb32260fd27ebbc1))

* Cleaned-up and updated some docstrings. ([`992eefb`](https://github.com/ollinevalainen/satellitetools/commit/992eefba8150073fbc69ba1ce30782a6fe4aaea9))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`49e6c46`](https://github.com/ollinevalainen/satellitetools/commit/49e6c468aa03cbdf0f5f709a91e16b300ea75a19))

* Update README.md ([`1e7e427`](https://github.com/ollinevalainen/satellitetools/commit/1e7e4277b9aa58052de7c43d77c444f59d8b66ff))

* Update README.md ([`672bde1`](https://github.com/ollinevalainen/satellitetools/commit/672bde16d4a303e9d702bbe608705b6d2ca87274))

* Update README.md ([`426780b`](https://github.com/ollinevalainen/satellitetools/commit/426780b4947a98976685e567de350bda042d7de1))

* Merge branch 'develop' of github.com:ollinevalainen/satellitetools into develop ([`21ff39c`](https://github.com/ollinevalainen/satellitetools/commit/21ff39c9e7adba7462e12b1b0fabfada25ce5c4f))

* Fixed bugs ([`1c1e5c0`](https://github.com/ollinevalainen/satellitetools/commit/1c1e5c05de680e9cc70bc71bf916f6314184cb01))

* Updated example ([`6e3ca22`](https://github.com/ollinevalainen/satellitetools/commit/6e3ca22a48f6966da6bc0d3769b0044f22574ade))

* Cleaning-up ([`da2af52`](https://github.com/ollinevalainen/satellitetools/commit/da2af529d8e9b620b6d3cccf90262da982f94e28))

* Switched xarray order to time, band, y, x ([`3c5e265`](https://github.com/ollinevalainen/satellitetools/commit/3c5e26524814ae08456c521b819d80245790e5ad))

* Removed added flipping of y-coords ([`2796497`](https://github.com/ollinevalainen/satellitetools/commit/2796497c96b779238e89be93d669c73e8c5579d9))

* Updated gee, aws_cog functions to return properly ([`2caedd2`](https://github.com/ollinevalainen/satellitetools/commit/2caedd2d5eb4dccc31de0a55d3d0c8136dac2656))

* Restrucured and created aws_cog pipeline ([`9588df7`](https://github.com/ollinevalainen/satellitetools/commit/9588df76c39d8006a70e4a154b1ed503d114e0bd))

* Update README.md ([`21010fe`](https://github.com/ollinevalainen/satellitetools/commit/21010fe0e8fe9da1e007e9ffeab9c04ff7e8639f))

* Autoformated with black ([`0c7f214`](https://github.com/ollinevalainen/satellitetools/commit/0c7f214f99ed4911ec6a43c596b4509364f3696e))

* Fixed tile handling/filtering ([`b37bf18`](https://github.com/ollinevalainen/satellitetools/commit/b37bf184c2512b508c2e9f44b56864150540d8cf))

* Added gitignore ([`5fc5f4d`](https://github.com/ollinevalainen/satellitetools/commit/5fc5f4d2857452979f50340ac6a7044b051aebf9))

* Added comments. ([`428cb01`](https://github.com/ollinevalainen/satellitetools/commit/428cb01b6c9d9e36b7cdfc405eb34146851d7f7b))

* Added kml/kmz reader. ([`9ea4e98`](https://github.com/ollinevalainen/satellitetools/commit/9ea4e98e589381311c00276c5b6334921c6ca79b))

* Added comments and cleaning-up. ([`a8a2d4b`](https://github.com/ollinevalainen/satellitetools/commit/a8a2d4bb9a0fb18cd02618318018d794eb522a7a))

* Merge pull request #1 from ollinevalainen/add-license-1

Create LICENSE ([`9cb781e`](https://github.com/ollinevalainen/satellitetools/commit/9cb781e12a70f7c4b62903d52afdf609fb6645fe))

* Create LICENSE ([`5ef6cfa`](https://github.com/ollinevalainen/satellitetools/commit/5ef6cfa0cd8f4344f265579374a03af3945fec41))

* Fixed paths. ([`dfefdfa`](https://github.com/ollinevalainen/satellitetools/commit/dfefdfa2cfd19fa46a2bf075b046d064e6a9ce35))

* Added output dir ([`2a7d52f`](https://github.com/ollinevalainen/satellitetools/commit/2a7d52f358ee6f4cf71d4481d80a7cba4d42ea45))

* Resctructured repo ([`04f043e`](https://github.com/ollinevalainen/satellitetools/commit/04f043eda4820a49fe0fe42e9bbd66178c9e81bb))

* Update README.md ([`bd34c34`](https://github.com/ollinevalainen/satellitetools/commit/bd34c34cb93b22e9593cbae35d5f81dc30c9e10a))

* Removed old file. ([`cd97639`](https://github.com/ollinevalainen/satellitetools/commit/cd976392c8d9d0e1bc63edf8889206a39dac3f13))

* Create README.md ([`a018a05`](https://github.com/ollinevalainen/satellitetools/commit/a018a05ba3604160f4f6a03638e145cce6fa1162))

* First commit. ([`392924a`](https://github.com/ollinevalainen/satellitetools/commit/392924aaa39eaa8aea6cb4f78f1928f26f677f1c))
